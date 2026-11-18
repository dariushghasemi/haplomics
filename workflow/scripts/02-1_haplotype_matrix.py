#!/usr/bin/env python3
"""Construct diploid haplotypes from ordered CHROM POS REF ALT GT... records.

Requires phased, non-missing genotypes at every locus variant for a sample.
For pipeline testing only, --allow-unphased-test sorts '/' alleles at
each site to construct synthetic haplotypes; it does not infer phase.
Frequency calculation precedes optional rare-carrier dosage filtering.
Uses Python standard library; genotype rows are streamed to byte buffers.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import platform
import sys


def checksum(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def construct(args):
    allow_unphased = getattr(args, "allow_unphased_test", False)
    construction_mode = "unphased_test" if allow_unphased else "phased"
    if not 0 <= args.min_freq <= 1:
        raise ValueError("--min-freq must be between zero and one")
    samples = [line.strip() for line in args.samples.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
    if not samples or len(set(samples)) != len(samples):
        raise ValueError("Sample IDs must be non-empty and unique")
    first, second = [bytearray() for _ in samples], [bytearray() for _ in samples]
    exclusion = [None] * len(samples)
    unphased_calls = [0] * len(samples)
    variants, seen = [], set()
    with args.gt.open(encoding="utf-8-sig", newline="") as handle:
        for line_number, fields in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not fields:
                continue
            if len(fields) != 4+len(samples):
                raise ValueError(f"GT line {line_number}: expected CHROM POS REF ALT plus {len(samples)} GT columns; got {len(fields)} columns")
            chrom, pos, ref, alt = fields[:4]
            if not all((chrom, ref, alt)) or not pos.isdigit() or int(pos) < 1:
                raise ValueError(f"Malformed variant on GT line {line_number}")
            if ref == "." or alt == ".":
                raise ValueError(f"Missing REF/ALT on GT line {line_number}")
            identifier = f"{chrom}:{pos}:{ref}:{alt}"
            if identifier in seen:
                raise ValueError(f"Duplicate variant: {identifier}")
            seen.add(identifier)
            variants.append((len(variants)+1,chrom,int(pos),ref,alt,identifier))
            allele_max = len(alt.split(","))
            for index, gt in enumerate(fields[4:]):
                if "." in gt:
                    reason = "missing_gt"
                elif "/" in gt and not allow_unphased:
                    reason = "unphased_gt"
                else:
                    if "/" in gt and "|" in gt:
                        raise ValueError(f"Invalid GT {gt!r} for {samples[index]} at {identifier}")
                    alleles = gt.replace("/", "|").split("|")
                    reason = "non_diploid_gt" if len(alleles) != 2 else None
                    if reason is None:
                        if any(not a.isdigit() or int(a)>allele_max for a in alleles):
                            raise ValueError(f"Invalid GT {gt!r} for {samples[index]} at {identifier}")
                        if "/" in gt:
                            # Unphased allele order has no chromosome meaning.
                            # Canonicalize numerically for a reproducible test fixture.
                            alleles.sort(key=int)
                            unphased_calls[index] += 1
                if reason:
                    if exclusion[index] is None:
                        exclusion[index] = reason
                        first[index].clear()
                        second[index].clear()
                    continue
                if exclusion[index] is None:
                    if len(variants)>1:
                        first[index].extend(b"-")
                        second[index].extend(b"-")
                    first[index].extend(str(int(alleles[0])).encode("ascii"))
                    second[index].extend(str(int(alleles[1])).encode("ascii"))
    if not variants:
        raise ValueError("No variants found")
    pairs, counts = [], Counter()
    for index in range(len(samples)):
        pair = None if exclusion[index] else (first[index].decode("ascii"),second[index].decode("ascii"))
        pairs.append(pair)
        if pair:
            counts.update(pair)
        first[index] = second[index] = None
    total = sum(counts.values())
    if not total:
        if allow_unphased:
            raise ValueError("No usable diploid genotypes; expected non-missing GT values such as 0|1 or 0/1")
        raise ValueError("No usable phased diploid haplotypes; expected GT values such as 0|1. "
                         "For unphased pipeline tests only, use --allow-unphased-test")
    if allow_unphased:
        print("TEST ONLY: '/' alleles are sorted per site to construct synthetic haplotypes; "
              "phase is not inferred. Do not use these results for haplotype association analysis.",file=sys.stderr)
    ordered = sorted(counts,key=lambda h:(-counts[h],h))
    # Assign once, before filtering: every output uses the same locus-local IDs.
    haplotype_ids = {h: f"H{i}" for i, h in enumerate(ordered, 1)}
    common = [h for h in ordered if counts[h]/total >= args.min_freq]
    if not common:
        raise ValueError("No haplotypes passed --min-freq")
    common_set = set(common)
    reference, tested = common[0], common[1:]
    prefix = str(args.out_prefix)
    args.out_prefix.parent.mkdir(parents=True,exist_ok=True)

    def write(suffix, header, records):
        with open(prefix+suffix,"w",encoding="utf-8",newline="") as handle:
            writer = csv.writer(handle,delimiter="\t",lineterminator="\n")
            writer.writerow(header)
            writer.writerows(records)

    write(".variants.tsv",["variant_index","chrom","pos","ref","alt","variant_id"],variants)
    # This complete table is also the ID -> ordered allele composition dictionary.
    # Rare compositions keep distinct IDs; only their analysis group is pooled.
    write(".frequencies.tsv",["haplotype","composition","count","frequency","status","group","construction_mode"],
          ((haplotype_ids[h],h,counts[h],f"{counts[h]/total:.12g}",
            "reference" if h==reference else "common_tested" if h in common_set else "rare",
            haplotype_ids[h] if h in common_set else "Rare",construction_mode) for h in ordered))
    write(".sample_haplotypes.tsv",["sample_id","hap1","hap2","usable"],
          ((s,*(tuple(haplotype_ids[h] for h in pair) if pair else ("NA","NA")),int(pair is not None))
           for s,pair in zip(samples,pairs)))
    dosage_columns = [haplotype_ids[h] for h in tested]
    if not args.drop_rare_carriers:
        dosage_columns.append("Rare")
    def dosage_rows():
        for sample,pair in zip(samples,pairs):
            if pair is None:
                yield [sample]+["NA"]*len(dosage_columns)
                continue
            rare_count = sum(h not in common_set for h in pair)
            if args.drop_rare_carriers and rare_count:
                continue
            yield [sample]+[pair.count(h) for h in tested]+([] if args.drop_rare_carriers else [rare_count])
    write(".dosage.tsv",["sample_id"]+dosage_columns,dosage_rows())
    metadata = dict(schema_version=3,
                    method="synthetic_unphased_test" if allow_unphased else "phased_diploid_extraction",
                    construction_mode=construction_mode,test_only=allow_unphased,
                    unphased_handling="sort allele indices numerically per '/' GT; preserve '|' order" if allow_unphased else "exclude sample",
                    n_unphased_genotypes_in_usable_samples=sum(n for n,pair in zip(unphased_calls,pairs) if pair is not None),
                    python_version=platform.python_version(),
                    inputs={str(p):checksum(p) for p in (args.gt,args.samples)},
                    script_sha256=checksum(__file__),variant_order="input_record_order",variant_index_base=1,
                    n_variants=len(variants),n_samples=len(samples),n_usable_samples=total//2,
                    n_excluded_samples=sum(x is not None for x in exclusion),
                    exclusions_by_first_reason=dict(Counter(x for x in exclusion if x)),
                    frequency_denominator=total,
                    frequency_denominator_definition=("2 * samples with non-missing, diploid GT at every variant; synthetic test copies"
                        if allow_unphased else "2 * samples with phased, non-missing, diploid GT at every variant"),
                    haplotype_naming="H1, H2, ... by descending count, then lexicographic composition",
                    haplotype_id_scope="local to this locus and dataset; assigned before frequency filtering",
                    haplotype_dictionary=Path(prefix+".frequencies.tsv").name,
                    rare_group="Rare",
                    min_frequency=args.min_freq,reference_haplotype=haplotype_ids[reference],n_common_haplotypes=len(common),
                    n_rare_haplotypes=len(counts)-len(common),drop_rare_carriers=args.drop_rare_carriers,
                    n_dosage_samples=sum(pair is None or not args.drop_rare_carriers or all(h in common_set for h in pair) for pair in pairs))
    Path(prefix+".metadata.json").write_text(json.dumps(metadata,indent=2)+"\n",encoding="utf-8")
    print(f"Reference: {haplotype_ids[reference]}; {len(variants)} variants; {total//2}/{len(samples)} usable samples",file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gt",required=True,type=Path,help="Headerless CHROM POS REF ALT GT... TSV")
    parser.add_argument("--samples",required=True,type=Path)
    parser.add_argument("--min-freq",type=float,default=0.01)
    parser.add_argument("--out-prefix",type=Path,default=Path("haplotypes"))
    parser.add_argument("--drop-rare-carriers",action="store_true")
    parser.add_argument("--allow-unphased-test",action="store_true",
                        help="Accept '/' GTs by sorting alleles per site; creates synthetic test haplotypes, not inferred phase")
    args = parser.parse_args()
    try:
        construct(args)
    except (ValueError,OSError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
