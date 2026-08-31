#!/usr/bin/env python3
"""Sequence-based filtering of ProteinMPNN designs.

Reads ProteinMPNN FASTA (or a plain one-sequence-per-line text file), removes
duplicates, then applies an ordered list of filters defined in a JSON config.

Positions are specified by *Nrf2 residue number*, not by column index, so a
single config works for every peptide length (10-16 aa). The peptide is
anchored at its C-terminus (Leu84 by default):

    column (1-based) = length - (anchor_residue - residue_number)

Usage:
    python filter_sequences.py --input 4IFL_10.fa --config ../inputs/keap1_filters.json \
                               --output 4IFL_10_filtered.txt --summary counts.csv
    python filter_sequences.py --self-test
"""

import argparse
import csv
import json
import sys
from pathlib import Path


def read_sequences(path):
    """Read a ProteinMPNN .fa (line after each '>') or a plain sequence list."""
    lines = Path(path).read_text().splitlines()
    if any(line.startswith(">") for line in lines):
        # The first record is ProteinMPNN's native/template sequence; it is kept.
        return [lines[i + 1].strip() for i, line in enumerate(lines)
                if line.startswith(">") and i + 1 < len(lines)]
    return [line.strip() for line in lines if line.strip()]


def column_of(residue_number, length, anchor_residue):
    """1-based column of a residue number in a C-terminally anchored peptide."""
    return length - (anchor_residue - residue_number)


def apply_step(seqs, step, anchor_residue):
    """Return the sequences surviving one filter step."""
    kind = step["type"]

    if kind == "motif":
        motifs = step["deny"]
        return [s for s in seqs if not any(m in s for m in motifs)]

    if kind == "position":
        allow = set(step["allow"]) if "allow" in step else None
        deny = set(step["deny"]) if "deny" in step else None
        if (allow is None) == (deny is None):
            raise ValueError(f"step {step['name']!r}: give exactly one of allow/deny")

        def ok(seq):
            for res in step["residues"]:
                col = column_of(res, len(seq), anchor_residue)
                if not 1 <= col <= len(seq):
                    return False        # residue not present at this peptide length
                aa = seq[col - 1]
                if allow is not None and aa not in allow:
                    return False
                if deny is not None and aa in deny:
                    return False
            return True

        return [s for s in seqs if ok(s)]

    raise ValueError(f"unknown step type: {kind}")


def run_filters(seqs, config, verbose=True):
    """Apply every step in order. Returns (survivors, [(step_name, count), ...])."""
    anchor = config["anchor_residue"]
    counts = [("num_seq", len(seqs))]

    seqs = sorted(set(seqs))          # sorted -> deterministic output across runs
    counts.append(("unique", len(seqs)))
    if verbose:
        print(f"{'unique':<28} {len(seqs):>8}")

    for step in config["steps"]:
        before = len(seqs)
        seqs = apply_step(seqs, step, anchor)
        counts.append((step["name"], len(seqs)))
        if verbose:
            print(f"{step['name']:<28} {len(seqs):>8}  (-{before - len(seqs)})")

    return seqs, counts


def self_test():
    """Minimal check of the residue->column mapping and both filter types."""
    # 15mer FFAQLQLDEETGEFL: ETGE core sits at columns 10-13 (Nrf2 E79-T80-G81-E82)
    assert column_of(79, 15, 84) == 10
    assert column_of(82, 15, 84) == 13
    assert column_of(79, 16, 84) == 11      # same residue, one column later in a 16mer
    assert column_of(84, 10, 84) == 10      # anchor is always the last column

    cfg = {"anchor_residue": 84, "steps": [
        {"name": "incl G81", "type": "position", "residues": [81], "allow": "G"},
        {"name": "excl ETGE", "type": "motif", "deny": ["ETGE"]},
    ]}
    seqs = ["FFAQLQLDEETGEFL",   # G at 81, but contains ETGE -> dropped by step 2
            "FFAQLQLDEEYGEFL",   # G at 81 (col 12), no ETGE          -> kept
            "FFAQLQLDEEYAEFL"]   # A at 81                            -> dropped by step 1
    kept, counts = run_filters(seqs, cfg, verbose=False)
    assert kept == ["FFAQLQLDEEYGEFL"], kept
    assert counts == [("num_seq", 3), ("unique", 3), ("incl G81", 2), ("excl ETGE", 1)], counts
    print("self-test OK")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, help="ProteinMPNN .fa or plain sequence list")
    ap.add_argument("--config", type=Path, help="JSON filter definition")
    ap.add_argument("--output", type=Path, help="surviving sequences, one per line")
    ap.add_argument("--summary", type=Path, help="optional CSV of per-step counts")
    ap.add_argument("--self-test", action="store_true", help="run internal checks and exit")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return

    if not (args.input and args.config and args.output):
        ap.error("--input, --config and --output are required")

    config = json.loads(args.config.read_text())
    seqs, counts = run_filters(read_sequences(args.input), config)

    args.output.write_text("\n".join(seqs) + "\n")
    print(f"\nwrote {len(seqs)} sequences -> {args.output}")

    if args.summary:
        with open(args.summary, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["step", "remaining"])
            w.writerows(counts)
        print(f"wrote step counts -> {args.summary}")


if __name__ == "__main__":
    sys.exit(main())
