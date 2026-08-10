#!/usr/bin/env python3
"""Generate child-reference informative-site tables for candidate trio NAHR events."""

from __future__ import annotations

import argparse
import re
import subprocess
from collections import OrderedDict
from pathlib import Path

GENE_UNITS = {"2", "3", "4", "6", "7"}
INTERGENIC_UNIT = "5"
VALID_BASES = {"A", "C", "G", "T"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bed", required=True, type=Path, help="Curated principal-bundle BED.")
    p.add_argument("--pairs", required=True, type=Path, help="TSV with child_hap and parent_hap columns.")
    p.add_argument("--fasta-root", required=True, type=Path, help="Root containing haplotype FASTA files.")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--mafft", default="mafft")
    p.add_argument("--threads", default="-1")
    p.add_argument("--snp-only", action="store_true")
    return p.parse_args()


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._#=-]+", "_", text)


def read_pairs(path: Path):
    import pandas as pd

    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"child_hap", "parent_hap"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")
    return list(df[["child_hap", "parent_hap"]].itertuples(index=False, name=None))


def read_blocks(path: Path, haplotype: str):
    records = []
    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"{path}:{line_no} has fewer than 6 columns")
            if fields[4] != haplotype:
                continue
            records.append((int(fields[1]), int(fields[2]), fields[5]))

    records.sort()
    blocks = []
    gene_no = 1
    inter_no = 1
    i = 0
    while i < len(records):
        start, end, unit = records[i]
        if unit in GENE_UNITS:
            j = i + 1
            while j < len(records) and records[j][2] in GENE_UNITS:
                j += 1
            blocks.append(
                {
                    "name": f"{haplotype}_Gene_{gene_no}",
                    "type": "Gene",
                    "start": min(x[0] for x in records[i:j]),
                    "end": max(x[1] for x in records[i:j]),
                }
            )
            gene_no += 1
            i = j
        elif unit == INTERGENIC_UNIT:
            blocks.append(
                {
                    "name": f"{haplotype}_Inter_{inter_no}",
                    "type": "Intergenic",
                    "start": start,
                    "end": end,
                }
            )
            inter_no += 1
            i += 1
        else:
            i += 1
    return blocks


def resolve_fasta(root: Path, haplotype: str) -> Path:
    sample = haplotype.split(".")[0]
    candidates = [
        root / sample / f"{haplotype}.PGA.fa",
        root / sample / f"{haplotype}.fa",
        root / f"{haplotype}.PGA.fa",
        root / f"{haplotype}.fa",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f"No FASTA found for {haplotype}; checked: {candidates}")


def read_single_fasta(path: Path) -> str:
    seq = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                seq.append(line.strip())
    if not seq:
        raise ValueError(f"No sequence found in {path}")
    return "".join(seq).upper()


def extract_blocks(blocks, sequence: str, block_type: str):
    seqs = OrderedDict()
    for block in blocks:
        if block["type"] != block_type:
            continue
        seq = sequence[block["start"] : block["end"]]
        if seq:
            seqs[block["name"]] = seq
    return seqs


def write_fasta(records: OrderedDict[str, str], path: Path) -> None:
    with path.open("w") as out:
        for name, seq in records.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")


def read_alignment(path: Path) -> OrderedDict[str, str]:
    records: OrderedDict[str, str] = OrderedDict()
    current = None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = ""
            else:
                if current is None:
                    raise ValueError(f"Invalid FASTA: {path}")
                records[current] += line.upper()
    lengths = {len(seq) for seq in records.values()}
    if len(lengths) != 1:
        raise ValueError(f"Alignment sequences have unequal lengths: {path}")
    return records


def run_mafft(records: OrderedDict[str, str], prefix: Path, mafft: str, threads: str):
    input_fa = prefix.with_suffix(".mafft_input.fa")
    output_fa = prefix.with_suffix(".mafft.fa")
    write_fasta(records, input_fa)
    cmd = [mafft, "--auto", "--thread", str(threads), "--quiet", str(input_fa)]
    with output_fa.open("w") as out:
        subprocess.run(cmd, check=True, stdout=out)
    input_fa.unlink()
    return read_alignment(output_fa), output_fa


def position_map(ref: str):
    positions = []
    pos = 0
    for base in ref:
        if base != "-":
            pos += 1
        positions.append(pos)
    return positions


def write_informative_vcf(
    alignment: OrderedDict[str, str],
    parent_blocks: list[str],
    ref_child: str,
    block_type: str,
    output: Path,
    snp_only: bool,
) -> int:
    samples = [name for name in alignment if name != ref_child]
    ref = alignment[ref_child]
    pos_map = position_map(ref)
    rows = []

    for i, ref_base in enumerate(ref):
        if ref_base == "-":
            continue
        if snp_only and ref_base not in VALID_BASES:
            continue

        parent_bases = {name: alignment[name][i] for name in parent_blocks}
        if any(base in {"-", "N"} for base in parent_bases.values()):
            continue
        if snp_only and any(base not in VALID_BASES for base in parent_bases.values()):
            continue
        if len(set(parent_bases.values())) <= 1:
            continue
        if ref_base not in set(parent_bases.values()):
            continue

        sample_bases = [alignment[name][i] for name in samples]
        if snp_only and any(base not in VALID_BASES for base in sample_bases):
            continue
        alts = sorted(set(sample_bases) - {ref_base})
        if not alts:
            continue

        allele_index = {ref_base: "0", **{alt: str(j + 1) for j, alt in enumerate(alts)}}
        genotypes = [allele_index.get(base, ".") for base in sample_bases]
        matched = [name for name, base in parent_bases.items() if base == ref_base]
        info = ";".join(
            [
                f"ALN_POS={i + 1}",
                f"BLOCK_TYPE={block_type}",
                f"REF_CHILD_BLOCK={ref_child}",
                f"MATCH_PARENT_BLOCKS={','.join(matched)}",
                "PARENT_BASES=" + ",".join(f"{name}:{base}" for name, base in parent_bases.items()),
            ]
        )
        rows.append(
            [ref_child, str(pos_map[i]), ".", ref_base, ",".join(alts), ".", "PASS", info, "GT", *genotypes]
        )

    with output.open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=nahr_mafft_to_childref_vcf.py\n")
        out.write("##note=Local pseudo-VCF for manual block-level NAHR inspection.\n")
        out.write("##coordinate=POS is the 1-based ungapped coordinate of the child reference block.\n")
        header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
    return len(rows)


def analyze_pair(child_hap: str, parent_hap: str, args: argparse.Namespace):
    child_blocks = read_blocks(args.bed, child_hap)
    parent_blocks = read_blocks(args.bed, parent_hap)
    child_seq = read_single_fasta(resolve_fasta(args.fasta_root, child_hap))
    parent_seq = read_single_fasta(resolve_fasta(args.fasta_root, parent_hap))

    pair_dir = args.outdir / f"{safe_name(child_hap)}__{safe_name(parent_hap)}"
    pair_dir.mkdir(parents=True, exist_ok=True)
    summary = []

    for block_type in ("Gene", "Intergenic"):
        parent_records = extract_blocks(parent_blocks, parent_seq, block_type)
        child_records = extract_blocks(child_blocks, child_seq, block_type)
        if len(parent_records) < 2 or not child_records:
            continue

        records = OrderedDict()
        records.update(parent_records)
        records.update(child_records)
        prefix = pair_dir / f"{safe_name(child_hap)}__{safe_name(parent_hap)}__{block_type}"
        alignment, alignment_path = run_mafft(records, prefix, args.mafft, args.threads)

        for child_block in child_records:
            output = pair_dir / f"{prefix.name}__ref-{safe_name(child_block)}.informative.vcf"
            n_sites = write_informative_vcf(
                alignment,
                list(parent_records),
                child_block,
                block_type,
                output,
                args.snp_only,
            )
            if n_sites == 0:
                output.unlink()
            summary.append(
                {
                    "child_hap": child_hap,
                    "parent_hap": parent_hap,
                    "block_type": block_type,
                    "ref_child_block": child_block,
                    "n_parent_blocks": len(parent_records),
                    "n_child_blocks": len(child_records),
                    "n_informative_sites": n_sites,
                    "alignment": str(alignment_path),
                    "vcf_like": str(output) if n_sites else ".",
                }
            )
    return summary


def main() -> None:
    import pandas as pd

    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    pairs = read_pairs(args.pairs)
    rows = []
    for child_hap, parent_hap in pairs:
        rows.extend(analyze_pair(child_hap, parent_hap, args))
    pd.DataFrame(rows).to_csv(args.outdir / "run_summary.tsv", sep="\t", index=False)
    print(f"Candidate pairs: {len(pairs)}")
    print(f"Wrote: {args.outdir / 'run_summary.tsv'}")


if __name__ == "__main__":
    main()
