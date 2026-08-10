#!/usr/bin/env python3

import argparse
import sys
from collections import Counter

from Bio import AlignIO


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize variable alignment sites among classified PGA paralog groups. "
            "Sequence IDs must begin with '<group>#'."
        )
    )
    parser.add_argument("--alignment", required=True, help="Input FASTA multiple-sequence alignment.")
    parser.add_argument("--output", required=True, help="Output VCF-like variable-site table.")
    parser.add_argument(
        "--reference",
        required=True,
        help="Sequence ID used as the alignment-coordinate reference.",
    )
    parser.add_argument(
        "--groups",
        nargs="+",
        required=True,
        help="Paralog groups to compare, e.g. PGA34A PGA34B PGA5.",
    )
    return parser.parse_args()


def group_from_id(sequence_id):
    return sequence_id.split("#", 1)[0]


def choose_major_alt(column, ref):
    """Return the most abundant non-reference allele; break ties lexicographically."""
    counts = Counter(base for base in column if base != ref)
    if not counts:
        return None
    return sorted(counts, key=lambda base: (-counts[base], base))[0]


def summarize_group(column, ref, alt):
    """Return the predominant REF/ALT state and its within-group frequency."""
    ref_count = column.count(ref)
    alt_count = column.count(alt)

    # Preserve the historical analysis behavior: ties are assigned to ALT.
    if alt_count >= ref_count:
        return "1", alt_count / len(column)
    return "0", ref_count / len(column)


def format_specific(groups, group_gt):
    states = sorted(set(group_gt.values()))

    if len(states) == 1:
        return -1, "-1"

    groups_by_state = {
        state: [group for group in groups if group_gt[group] == state]
        for state in states
    }

    # GT is binary here, so at most two predominant states are expected.
    if len(states) != 2:
        raise ValueError(f"Unexpected predominant states: {states}")

    left, right = states
    if len(groups_by_state[left]) == len(groups_by_state[right]):
        return (
            len(groups_by_state[left]),
            f"{','.join(groups_by_state[left])}|{','.join(groups_by_state[right])}",
        )

    minority_state = min(states, key=lambda state: len(groups_by_state[state]))
    minority_groups = groups_by_state[minority_state]
    return len(minority_groups), ",".join(minority_groups)


def main():
    args = parse_args()

    if len(set(args.groups)) != len(args.groups):
        raise ValueError("--groups contains duplicate labels.")

    alignment = AlignIO.read(args.alignment, "fasta")
    records = list(alignment)
    sequence_ids = [record.id for record in records]

    if len(sequence_ids) != len(set(sequence_ids)):
        raise ValueError("The alignment contains duplicate sequence IDs.")
    if args.reference not in sequence_ids:
        raise ValueError(f"Reference sequence '{args.reference}' was not found.")

    group_indices = {}
    for group in args.groups:
        indices = [
            i
            for i, sequence_id in enumerate(sequence_ids)
            if group_from_id(sequence_id) == group
        ]
        if not indices:
            raise ValueError(
                f"No sequences found for group '{group}'. "
                "Expected IDs to begin with '<group>#'."
            )
        group_indices[group] = indices

    focal_indices = sorted(
        {index for indices in group_indices.values() for index in indices}
    )

    reference_index = sequence_ids.index(args.reference)
    reference_sequence = str(records[reference_index].seq).upper()
    alignment_length = len(reference_sequence)

    rows = []

    for i in range(alignment_length):
        column = [str(records[index].seq[i]).upper() for index in focal_indices]

        if "N" in column:
            print(
                f"WARNING: requested groups contain N at alignment position {i + 1}",
                file=sys.stderr,
            )

        if len(set(column)) <= 1:
            continue

        ref = reference_sequence[i]
        alt = choose_major_alt(column, ref)
        if alt is None:
            continue

        group_gt = {}
        group_values = {}

        for group in args.groups:
            group_column = [
                str(records[index].seq[i]).upper()
                for index in group_indices[group]
            ]
            gt, gf = summarize_group(group_column, ref, alt)
            group_gt[group] = gt
            group_values[group] = f"{gt}:{gf:.4f}"

        counts = Counter(column)
        ac = counts[alt]
        ns = len(column)
        af = ac / ns

        specific_count, specific = format_specific(args.groups, group_gt)
        info = (
            f"AC={ac};AF={af:.4f};NS={ns};"
            f"SPECIFIC_COUNT={specific_count};SPECIFIC={specific}"
        )

        rows.append(
            [
                args.reference,
                str(i + 1),
                ".",
                ref,
                alt,
                ".",
                "PASS",
                info,
                "GT:GF",
                *[group_values[group] for group in args.groups],
            ]
        )

    with open(args.output, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="ALT allele count across requested groups">\n'
        )
        out.write(
            '##INFO=<ID=AF,Number=A,Type=Float,Description="ALT allele frequency across requested groups">\n'
        )
        out.write(
            '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of sequences across requested groups">\n'
        )
        out.write(
            '##INFO=<ID=SPECIFIC_COUNT,Number=1,Type=Integer,Description="Number of paralog groups carrying the minority REF/ALT state">\n'
        )
        out.write(
            '##INFO=<ID=SPECIFIC,Number=.,Type=String,Description="Paralog group(s) distinguished at this site">\n'
        )
        out.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Predominant REF/ALT state within the paralog group">\n'
        )
        out.write(
            '##FORMAT=<ID=GF,Number=1,Type=Float,Description="Frequency of the predominant REF/ALT state within the paralog group">\n'
        )
        out.write(f"##contig=<ID={args.reference},length={alignment_length}>\n")
        out.write(
            "\t".join(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                ]
                + args.groups
            )
            + "\n"
        )

        for row in rows:
            out.write("\t".join(row) + "\n")

    print(f"Alignment length: {alignment_length}", file=sys.stderr)
    print(f"Variable sites: {len(rows)}", file=sys.stderr)


if __name__ == "__main__":
    main()
