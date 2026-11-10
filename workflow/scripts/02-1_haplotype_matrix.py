#!/usr/bin/env python3

from collections import Counter
import argparse
import sys


# CHANGED 1: Removed safe_name(); short labels are assigned below.


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gt", required=True, help="phased_gt.tsv from bcftools query")
    parser.add_argument("--samples", required=True, help="samples.txt from bcftools query -l")
    parser.add_argument("--min-freq", type=float, default=0.01)
    parser.add_argument("--out-prefix", default="haplotypes")
    parser.add_argument(
        "--drop-rare-carriers",
        action="store_true",
        help="Drop individuals carrying any rare haplotype instead of using Rare dosage",
    )

    # CHANGED 2: Explicitly separate real phased input from unphased testing.
    parser.add_argument(
        "--allow-unphased-test",
        action="store_true",
        help="TEST ONLY: accept / genotypes in their written allele order; does not infer phase",
    )
    args = parser.parse_args()
    if not 0 <= args.min_freq <= 1:
        parser.error("--min-freq must be between 0 and 1")
    if args.allow_unphased_test:
        print(
            "WARNING: unphased test mode; written allele order is arbitrary. "
            "Outputs must not be interpreted as inferred haplotypes or their frequencies.",
            file=sys.stderr,
        )

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

            if len(fields) < 5:
                sys.exit(f"Malformed line: {line}")

            chrom = fields[0]
            pos = fields[1]
            ref = fields[2]
            alt = fields[3]
            gts = fields[4:]

            if len(gts) != len(samples):
                sys.exit(
                    f"Number of genotypes does not match number of samples at {chrom}:{pos}. "
                    f"Found {len(gts)} genotypes but {len(samples)} samples."
                )

            variants.append((chrom, pos, ref, alt))
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

            # CHANGED 3: Require | by default; permit / only in test mode.
            # Input must contain GT values only, not GT:DP:... fields.
            if "." in gt:
                usable = False
                break
            if "|" in gt and "/" not in gt:
                alleles = gt.split("|")
            elif args.allow_unphased_test and "/" in gt and "|" not in gt:
                alleles = gt.split("/")
            else:
                usable = False
                break
            if any(not a.isdigit() for a in alleles):
                usable = False
                break

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
        sys.exit("No usable haplotypes found. Default mode requires phased diploid GT values (|); use --allow-unphased-test only for output testing.")

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

    # CHANGED 4: Deterministic labels in decreasing count order.
    # The existing allele-pattern tie-break makes ties independent of sample order.
    # Labels are local to this dataset and ordered variant list, not universal IDs.
    common_set = set(common_haps)
    hap_labels = {hap: f"H{i}" for i, hap in enumerate(common_haps, start=1)}
    ordered_haps = sorted(hap_counter, key=lambda h: (-hap_counter[h], h))
    rare_count = sum(hap_counter[h] for h in ordered_haps if h not in common_set)

    def label(hap):
        return hap_labels.get(hap, "Rare")

    print(f"Reference haplotype: {label(reference_hap)} ({reference_hap})", file=sys.stderr)
    print(f"Number of common haplotypes: {len(common_haps)}", file=sys.stderr)
    print(f"Number of test haplotypes: {len(test_haps)}", file=sys.stderr)

    # Output variant list
    with open(args.out_prefix + ".variants.tsv", "w") as out:
        out.write("variant_index\tchrom\tpos\tref\talt\n")
        for i, (chrom, pos, ref, alt) in enumerate(variants, start=1):
            out.write(f"{i}\t{chrom}\t{pos}\t{ref}\t{alt}\n")

    # CHANGED 5: One summary row per common label and one pooled Rare row.
    # Frequencies retain the original denominator: all usable chromosomes,
    # including those from rare carriers when --drop-rare-carriers is selected.
    with open(args.out_prefix + ".frequencies.tsv", "w") as out:
        out.write("haplotype\tcomposition\tcount\tfrequency\tstatus\n")
        for hap in common_haps:
            status = "reference" if hap == reference_hap else "common_tested"
            out.write(
                f"{label(hap)}\t{hap}\t{hap_counter[hap]}\t{hap_freqs[hap]:.8f}\t{status}\n"
            )
        if rare_count:
            out.write(f"Rare\t\t{rare_count}\t{rare_count / total_haps:.8f}\trare\n")

    # CHANGED 6: Lossless lookup; Rare is a group, not a single allele pattern.
    # Pattern positions correspond exactly to .variants.tsv variant_index.
    with open(args.out_prefix + ".haplotype_key.tsv", "w") as out:
        out.write("haplotype\tallele_pattern\tcount\tfrequency\tstatus\n")
        for hap in ordered_haps:
            status = (
                "reference" if hap == reference_hap
                else "common_tested" if hap in common_set
                else "rare"
            )
            out.write(
                f"{label(hap)}\t{hap}\t{hap_counter[hap]}\t{hap_freqs[hap]:.8f}\t{status}\n"
            )

    # Output per-sample haplotype pair
    with open(args.out_prefix + ".sample_haplotypes.tsv", "w") as out:
        # CHANGED 7: Short labels plus original patterns preserve sample-level detail.
        out.write("sample_id\thap1\thap2\tusable\thap1_pattern\thap2_pattern\n")

        for sample in samples:
            h1, h2 = sample_haps[sample]

            if h1 is None:
                out.write(f"{sample}\tNA\tNA\t0\tNA\tNA\n")
            else:
                out.write(f"{sample}\t{label(h1)}\t{label(h2)}\t1\t{h1}\t{h2}\n")

    # Output dosage matrix
    # CHANGED 8: Dosage columns use H2, H3, ... and Rare; H1 stays omitted
    # as the reference category, just as in the original script.
    dosage_columns = [label(h) for h in test_haps]

    if not args.drop_rare_carriers:
        dosage_columns.append("Rare")

    with open(args.out_prefix + ".dosage.tsv", "w") as out:
        out.write("\t".join(["sample_id"] + dosage_columns) + "\n")

        for sample in samples:
            h1, h2 = sample_haps[sample]

            if h1 is None:
                out.write("\t".join([sample] + ["NA"] * len(dosage_columns)) + "\n")
                continue

            sample_hap_list = [h1, h2]

            carries_rare = any(h not in common_set for h in sample_hap_list)

            if args.drop_rare_carriers and carries_rare:
                continue

            values = []

            for hap in test_haps:
                values.append(str(sample_hap_list.count(hap)))

            if not args.drop_rare_carriers:
                rare_count = sum(h not in common_set for h in sample_hap_list)
                values.append(str(rare_count))

            out.write("\t".join([sample] + values) + "\n")


if __name__ == "__main__":
    main()
