#!/usr/bin/env python3
"""Plot constructed haplotypes using an ordered bcftools SNP list.

Requires Matplotlib. Writes PREFIX.pdf (all pages), PREFIX.png and PREFIX.svg
(first page), and PREFIX.manifest.json. Does not change analysis statuses.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
from pathlib import Path


def read_variant_identifiers(path, indexed=False):
    if not indexed:
        return [line.strip() for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
    identifiers = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"variant_index", "chrom", "pos", "ref", "alt"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Variant TSV requires variant_index, chrom, pos, ref, alt columns")
        for expected, row in enumerate(reader, 1):
            if int(row["variant_index"]) != expected:
                raise ValueError("variant_index must be consecutive, one-based, and in construction order")
            identifier = ":".join(row[k] for k in ("chrom", "pos", "ref", "alt"))
            if "variant_id" in row and row["variant_id"] != identifier:
                raise ValueError(f"variant_id does not match its allele columns at index {expected}")
            identifiers.append(identifier)
    return identifiers


def read_inputs(frequencies, snps, indexed_variants=False):
    variants, seen_variants = [], set()
    for line in read_variant_identifiers(snps, indexed=indexed_variants):
        if not line.strip():
            continue
        identifier = line.strip()
        fields = identifier.split(":", 3)
        if len(fields) != 4 or not all(fields):
            raise ValueError(f"Invalid SNP identifier: {identifier!r}; expected CHROM:POS:REF:ALT")
        try:
            if int(fields[1]) < 1:
                raise ValueError
        except ValueError:
            raise ValueError(f"Invalid SNP position: {identifier!r}") from None
        if identifier in seen_variants:
            raise ValueError(f"Duplicate SNP identifier: {identifier}")
        seen_variants.add(identifier)
        variants.append(dict(id=identifier, chrom=fields[0], pos=int(fields[1]),
                             ref=fields[2], alt=fields[3].split(",")))
    if not variants:
        raise ValueError("The SNP file is empty.")
    rows, seen, seen_compositions = [], set(), set()
    with frequencies.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"haplotype", "count", "frequency", "status"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Frequency TSV must contain: {', '.join(sorted(required))}")
        for line_number, record in enumerate(reader, 2):
            hap = record["haplotype"].strip()
            if hap in seen:
                raise ValueError(f"Duplicate haplotype on TSV line {line_number}: {hap}")
            seen.add(hap)
            # New tables carry short IDs and an explicit composition dictionary.
            # Older standalone tables put the composition in 'haplotype'.
            composition = record["composition"].strip() if "composition" in record else hap
            if "composition" in record and (not hap.startswith("H") or not hap[1:].isdigit()
                                             or int(hap[1:]) < 1 or hap != f"H{int(hap[1:])}"):
                raise ValueError(f"Invalid haplotype ID on TSV line {line_number}: {hap!r}")
            if composition in seen_compositions:
                raise ValueError(f"Duplicate composition on TSV line {line_number}: {composition}")
            seen_compositions.add(composition)
            tokens = composition.split("-")
            if len(tokens) != len(variants):
                raise ValueError(f"TSV line {line_number}: {len(tokens)} alleles but {len(variants)} SNPs")
            alleles = []
            for token, variant in zip(tokens, variants):
                if token == ".":
                    alleles.append(None)
                elif token.isdigit() and int(token) <= len(variant["alt"]):
                    alleles.append(int(token))
                else:
                    raise ValueError(f"Invalid allele {token!r} at {variant['id']} on TSV line {line_number}")
            count, frequency = float(record["count"]), float(record["frequency"])
            if not math.isfinite(count) or count < 0:
                raise ValueError(f"Invalid count on TSV line {line_number}")
            if not math.isfinite(frequency) or not 0 <= frequency <= 1:
                raise ValueError(f"Invalid frequency on TSV line {line_number}")
            status = record["status"].strip()
            construction_mode = record.get("construction_mode", "phased").strip()
            if construction_mode not in {"phased", "unphased_test"}:
                raise ValueError(f"Invalid construction_mode on TSV line {line_number}: {construction_mode!r}")
            identity = "\n".join(v["id"] for v in variants) + "\n" + composition
            rows.append(dict(haplotype=hap, composition=composition, alleles=alleles, count=count,
                             frequency=frequency, status=status,
                             construction_mode=construction_mode,
                             label=hap if "composition" in record else None,
                             id="H" + hashlib.sha256(identity.encode()).hexdigest()[:12]))
    if not rows or sum(r["count"] for r in rows) <= 0:
        raise ValueError("The frequency table must contain a positive total count.")
    if len({r["construction_mode"] for r in rows}) != 1:
        raise ValueError("Frequency table mixes phased and unphased_test construction modes")
    # Allow rounding of the frequency column (the example uses eight decimals).
    tolerance = max(1e-6, len(rows) * 5.1e-9)
    if abs(sum(r["frequency"] for r in rows) - 1) > tolerance:
        raise ValueError("Frequencies do not sum to one; provide the complete frequency table.")
    total = sum(r["count"] for r in rows)
    if any(abs(r["frequency"] - r["count"] / total) > 1e-6 for r in rows):
        raise ValueError("Counts and frequencies disagree (absolute tolerance 1e-6).")
    rows.sort(key=lambda r: (-r["count"], r["composition"]))
    for index, row in enumerate(rows, 1):
        if row["label"] is None:
            row["label"] = f"H{index}"
    return variants, rows


def select_rows(rows, maximum, minimum):
    eligible = [r for r in rows if r["status"] != "rare" and r["frequency"] >= minimum]
    chosen = eligible[:maximum] if maximum else eligible
    ids = {r["id"] for r in chosen}
    ids.update(r["id"] for r in rows if r["status"] in {"reference", "common_tested"})
    selected = [r for r in rows if r["id"] in ids]
    return sorted(selected, key=lambda r: r["status"] != "reference")


def count_label(value):
    return f"{value:,.0f}" if value.is_integer() else f"{value:,.3f}".rstrip("0").rstrip(".")


def allele_text(variant, allele):
    if allele is None:
        return "."
    return variant["ref"] if allele == 0 else variant["alt"][allele-1]


def display_alleles(variants, alleles, reference, is_reference=False):
    # Missing alleles remain explicit, including when the reference is missing.
    return [allele_text(v, a) if is_reference or a is None or r is None or a != r else ""
            for v, a, r in zip(variants, alleles, reference)]


def draw_page(plt, np, colors, variants, rows, all_rows, reference, args, page, pages,
              row_start, variant_start, shown_frequency):
    n, m = len(rows), len(variants)
    width = max(9, 4.5 + 0.50 * m)
    test_mode = args.construction_mode == "unphased_test"
    extra_header = 0.3 if test_mode else 0
    height = 3.3 + extra_header + 0.42 * n
    fig = plt.figure(figsize=(width, height), facecolor="white")
    grid = fig.add_gridspec(1, 2, width_ratios=[max(2.4, 0.50*m), 2.0],
                           left=0.25, right=0.97, bottom=1.05/height,
                           top=1-(2.05+extra_header)/height, wspace=0.12)
    ax = fig.add_subplot(grid[0, 0])
    table = fig.add_subplot(grid[0, 1], sharey=ax)
    ax.imshow(np.zeros((n,m)), cmap=colors.ListedColormap(["#edf1f5"]), aspect="auto",
              interpolation="nearest")
    ax.set_xticks(range(m), [v["id"] for v in variants], rotation=55, ha="left",
                  rotation_mode="anchor", fontsize=11)
    ax.xaxis.tick_top()
    labels = [r["label"] for r in rows]
    labels = [label + " (reference)" if r["status"] == "reference" else label
              for label,r in zip(labels,rows)]
    ax.set_yticks(range(n), labels, fontsize=11)
    ax.tick_params(axis="both", length=0, pad=5)
    ax.set_xticks(np.arange(-0.5, m), minor=True)
    ax.set_yticks(np.arange(-0.5, n), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.4)
    ax.tick_params(which="minor", length=0)
    reference_block = reference["alleles"][variant_start:variant_start+m]
    for i, row in enumerate(rows):
        texts = display_alleles(variants, row["alleles"], reference_block,
                               is_reference=row["id"] == reference["id"])
        for j, text in enumerate(texts):
            ax.text(j, i, text, ha="center", va="center", fontsize=13,
                    fontweight="bold" if row["id"] == reference["id"] else "normal",
                    color="#34495e")
    table.set_xlim(0, 1)
    table.axis("off")
    table.text(0.50, 1.05, "Count", transform=table.transAxes, ha="right", fontsize=12, fontweight="bold")
    table.text(0.98, 1.05, "%", transform=table.transAxes, ha="right", fontsize=12, fontweight="bold")
    for i, row in enumerate(rows):
        table.text(0.50, i, count_label(row["count"]), ha="right", va="center", fontsize=12)
        pct = row["frequency"]*100
        label = f"{pct:.3g}" if 0 < pct < 0.001 else f"{pct:.3f}"
        table.text(0.98, i, label, ha="right", va="center", fontsize=12)
    for axis in (ax, table):
        axis.set_ylim(n-0.5, -0.5)
        for spine in axis.spines.values():
            spine.set_visible(False)
    total = sum(r["count"] for r in all_rows)
    fig.suptitle(args.title or args.frequencies.name.removesuffix(".frequencies.tsv"),
                 x=0.03, y=1-0.15/height, ha="left", fontsize=16, fontweight="bold")
    fig.text(0.03, 1-0.51/height,
             f"{len(all_rows):,} distinct haplotypes | {count_label(total)} counted copies",
             fontsize=10, color="#526577")
    fig.text(0.03, 1-0.76/height,
             f"Page {page}/{pages}: rows {row_start+1}–{row_start+n}, variants {variant_start+1}–{variant_start+m}",
             fontsize=10, color="#526577")
    if test_mode:
        fig.text(0.03, 1-1.03/height,
                 "TEST ONLY: synthetic haplotypes from unphased GTs; phase was not inferred.",
                 fontsize=10, color="#a33b20", fontweight="bold")
    omitted = len(all_rows)-args.selected_count
    rare = [r for r in all_rows if r["status"] == "rare"]
    rare_pct = sum(r["frequency"] for r in rare)*100
    fig.text(0.03, 0.70/height,
             "Blank cell: same allele as the reference haplotype. Period (.): missing allele.",
             fontsize=10, color="#526577")
    fig.text(0.03, 0.44/height,
             f"Rare (pooled; omitted from allele rows): {len(rare):,} haplotypes, "
             f"{count_label(sum(r['count'] for r in rare))} copies ({rare_pct:.3f}%).",
             fontsize=10, color="#526577")
    fig.text(0.03, 0.18/height,
             f"Shown: {args.selected_count:,} ({shown_frequency*100:.3f}%). "
             f"Other omitted: {omitted-len(rare):,}. Percentages use all copies.",
             fontsize=10, color="#526577")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frequencies", required=True, type=Path)
    variant_source = parser.add_mutually_exclusive_group(required=True)
    variant_source.add_argument("--snps", type=Path,
                        help="Ordered, headerless CHROM:POS:REF:ALT list from bcftools query")
    variant_source.add_argument("--variants", type=Path,
                        help="Ordered .variants.tsv from the haplotype matrix script (preferred)")
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--title")
    parser.add_argument("--max-haplotypes", type=int, default=30,
                        help="Top N non-rare rows before adding reference/common_tested; 0 selects all non-rare rows")
    parser.add_argument("--min-frequency", type=float, default=0,
                        help="Minimum fraction for ordinary rows; reference/common_tested always included")
    parser.add_argument("--rows-per-page", type=int, default=30)
    parser.add_argument("--variants-per-page", type=int, default=30)
    parser.add_argument("--dpi", type=int, default=200)
    args = parser.parse_args()
    if (args.max_haplotypes < 0 or not 0 <= args.min_frequency <= 1 or
        args.rows_per_page < 1 or args.variants_per_page < 1 or args.dpi < 1):
        parser.error("Invalid selection, pagination, or resolution parameter")
    try:
        variant_path = args.variants if args.variants is not None else args.snps
        variants, rows = read_inputs(args.frequencies, variant_path, indexed_variants=args.variants is not None)
        args.construction_mode = rows[0]["construction_mode"]
        references = [r for r in rows if r["status"] == "reference"]
        if len(references) != 1:
            raise ValueError("Exactly one haplotype must have status 'reference' for allele comparison.")
        reference = references[0]
        selected = select_rows(rows, args.max_haplotypes, args.min_frequency)
        if not selected:
            raise ValueError("No haplotypes meet the selection criteria.")
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.error(str(error))
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    # Keep Matplotlib's cache inside the output folder on restricted workers.
    os.environ.setdefault("MPLCONFIGDIR", str(args.output_prefix.parent / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import colors
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams.update({"font.family":"DejaVu Sans", "svg.fonttype":"none", "pdf.fonttype":42,
                         "svg.hashsalt":"haplomics-haplotype-plot-v1"})
    args.selected_count = len(selected)
    shown_frequency = sum(r["frequency"] for r in selected)
    pages = math.ceil(len(selected)/args.rows_per_page)*math.ceil(len(variants)/args.variants_per_page)
    prefix = str(args.output_prefix)
    manifest = dict(schema_version=4, matplotlib_version=matplotlib.__version__,
                    construction_mode=args.construction_mode,test_only=args.construction_mode == "unphased_test",
                    python_version=platform.python_version(), numpy_version=np.__version__,
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    parameters={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
                    inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in (args.frequencies,variant_path)},
                    variants=variants, selected_haplotypes=selected,
                    omitted_haplotypes=len(rows)-len(selected),
                    dropped_rare_haplotypes=sum(r["status"] == "rare" for r in rows),
                    dropped_rare_frequency=sum(r["frequency"] for r in rows if r["status"] == "rare"),
                    reference_haplotype_id=reference["id"],
                    reference_haplotype_label=reference["label"],
                    selected_frequency=shown_frequency, preview_page=1, pages=[])
    with PdfPages(prefix+".pdf", metadata={"Title":args.title or args.frequencies.stem,
                                         "CreationDate":None,"ModDate":None}) as pdf:
        page = 0
        for row_start in range(0,len(selected),args.rows_per_page):
            group = selected[row_start:row_start+args.rows_per_page]
            for variant_start in range(0,len(variants),args.variants_per_page):
                page += 1
                block = variants[variant_start:variant_start+args.variants_per_page]
                sliced = [dict(r,alleles=r["alleles"][variant_start:variant_start+len(block)]) for r in group]
                fig = draw_page(plt,np,colors,block,sliced,rows,reference,args,page,pages,
                                row_start,variant_start,shown_frequency)
                pdf.savefig(fig)
                if page == 1:
                    fig.savefig(prefix+".png",dpi=args.dpi)
                    fig.savefig(prefix+".svg",metadata={"Date":None})
                plt.close(fig)
                manifest["pages"].append(dict(page=page, row_start=row_start+1,
                    row_end=row_start+len(group), variant_start=variant_start+1,
                    variant_end=variant_start+len(block)))
    Path(prefix+".manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print(f"Wrote {pages} PDF page(s), first-page PNG/SVG, and manifest to {prefix}.*")


if __name__ == "__main__":
    main()
