#!/usr/bin/env python3

from collections import Counter
import argparse
import re
import sys


def safe_name(hap):
    """
    Convert haplotype like 1-1-1 to column name H_1_1_1.
    """
    return "H_" + re.sub(r"[^A-Za-z0-9]+", "_", hap)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gt", required=True, help="phased_gt.tsv from bcftools query")
    parser.add_argument("--samples", required=True, help="samples.txt from bcftools query -l")
    parser.add_argument("--min-freq", type=float, default=0.01)
    parser.add_argument("--out-prefix", default="haplotypes")
    parser.add_argument(
        "--drop-rare-carriers",
        action="store_true",
        help="Drop individuals carrying any rare haplotype instead of using OTHER_RARE dosage",
    )

    args = parser.parse_args()

    # Read sample names
    with open(args.samples) as f:
        samples = [line.strip() for line in f if line.strip()]

    if not samples:
        sys.exit("No sample names found.")

    # Read genotype table
    rows = []
    variants = []

    with open(args.gt) as f:
        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 3:
                sys.exit(f"Malformed line: {line}")

            chrom = fields[0]
            pos = fields[1]
            gts = fields[2:]

            if len(gts) != len(samples):
                sys.exit(
                    f"Number of genotypes does not match number of samples at {chrom}:{pos}. "
                    f"Found {len(gts)} genotypes but {len(samples)} samples."
                )

            variants.append((chrom, pos))
            rows.append(gts)

    n_variants = len(rows)
    n_samples = len(samples)

    if n_variants == 0:
        sys.exit("No variants found.")

    # Build two haplotypes per sample
    sample_haps = {}
    hap_counter = Counter()

    for sample_idx, sample in enumerate(samples):
        hap1 = []
        hap2 = []
        usable = True

        for variant_idx in range(n_variants):
            gt = rows[variant_idx][sample_idx]

            # Require phased non-missing genotype
            if "/" not in gt or "." in gt:
                usable = False
                break

            alleles = gt.split("/")

            if len(alleles) != 2:
                usable = False
                break

            a1, a2 = alleles

            hap1.append(a1)
            hap2.append(a2)

        if usable:
            h1 = "-".join(hap1)
            h2 = "-".join(hap2)

            sample_haps[sample] = (h1, h2)

            hap_counter[h1] += 1
            hap_counter[h2] += 1
        else:
            sample_haps[sample] = (None, None)

    total_haps = sum(hap_counter.values())

    if total_haps == 0:
        sys.exit("No usable phased haplotypes found.")

    # Determine common haplotypes
    hap_freqs = {
        hap: count / total_haps
        for hap, count in hap_counter.items()
    }

    common_haps = [
        hap for hap, freq in hap_freqs.items()
        if freq >= args.min_freq
    ]

    common_haps = sorted(
        common_haps,
        key=lambda h: (-hap_counter[h], h)
    )

    if not common_haps:
        sys.exit("No haplotypes passed the frequency threshold.")

    reference_hap = common_haps[0]
    test_haps = common_haps[1:]

    print(f"Reference haplotype: {reference_hap}", file=sys.stderr)
    print(f"Number of common haplotypes: {len(common_haps)}", file=sys.stderr)
    print(f"Number of test haplotypes: {len(test_haps)}", file=sys.stderr)

    # Output variant list
    with open(args.out_prefix + ".variants.tsv", "w") as out:
        out.write("variant_index\tchrom\tpos\n")
        for i, (chrom, pos) in enumerate(variants, start=1):
            out.write(f"{i}\t{chrom}\t{pos}\n")

    # Output haplotype frequencies
    with open(args.out_prefix + ".frequencies.tsv", "w") as out:
        out.write("haplotype\tcount\tfrequency\tstatus\n")

        for hap, count in hap_counter.most_common():
            freq = count / total_haps

            if hap == reference_hap:
                status = "reference"
            elif hap in test_haps:
                status = "common_tested"
            else:
                status = "rare"

            out.write(f"{hap}\t{count}\t{freq:.8f}\t{status}\n")

    # Output per-sample haplotype pair
    with open(args.out_prefix + ".sample_haplotypes.tsv", "w") as out:
        out.write("sample_id\thap1\thap2\tusable\n")

        for sample in samples:
            h1, h2 = sample_haps[sample]

            if h1 is None:
                out.write(f"{sample}\tNA\tNA\t0\n")
            else:
                out.write(f"{sample}\t{h1}\t{h2}\t1\n")

    # Output dosage matrix
    dosage_columns = [safe_name(h) for h in test_haps]

    if not args.drop_rare_carriers:
        dosage_columns.append("OTHER_RARE")

    with open(args.out_prefix + ".dosage.tsv", "w") as out:
        out.write("sample_id\t" + "\t".join(dosage_columns) + "\n")

        for sample in samples:
            h1, h2 = sample_haps[sample]

            if h1 is None:
                out.write(sample + "\t" + "\t".join(["NA"] * len(dosage_columns)) + "\n")
                continue

            sample_hap_list = [h1, h2]

            carries_rare = any(h not in common_haps for h in sample_hap_list)

            if args.drop_rare_carriers and carries_rare:
                continue

            values = []

            for hap in test_haps:
                values.append(str(sample_hap_list.count(hap)))

            if not args.drop_rare_carriers:
                rare_count = sum(h not in common_haps for h in sample_hap_list)
                values.append(str(rare_count))

            out.write(sample + "\t" + "\t".join(values) + "\n")


if __name__ == "__main__":
    main()
