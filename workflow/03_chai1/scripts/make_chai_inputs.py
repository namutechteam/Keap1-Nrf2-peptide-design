#!/usr/bin/env python3
"""Split a filtered sequence list into per-sequence Chai-1 FASTA inputs.

Chai-1 folds one FASTA at a time, so each surviving peptide becomes its own file.
Files are dealt round-robin into N batch directories, one batch per GPU.

Usage:
    python make_chai_inputs.py --input 4IFL_10_filtered.txt --outdir chai_in --batches 2
    # -> chai_in/batch0/0.fasta, chai_in/batch1/1.fasta, ...
"""

import argparse
from pathlib import Path


def write_inputs(sequences, outdir, batches, receptor=None):
    """Write one FASTA per sequence, dealt round-robin across `batches` subdirs."""
    dirs = [outdir / f"batch{i}" for i in range(batches)]
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)

    for i, seq in enumerate(sequences):
        records = []
        if receptor:
            records += [">protein|name=receptor", receptor]
        records += [f">protein|name={i}", seq]
        (dirs[i % batches] / f"{i}.fasta").write_text("\n".join(records) + "\n")

    return len(sequences)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, required=True,
                    help="filtered sequences, one per line (output of step 02)")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--batches", type=int, default=1,
                    help="number of batch directories, typically one per GPU")
    ap.add_argument("--receptor-fasta", type=Path,
                    help="optional Keap1 Kelch-domain sequence, prepended to every input "
                         "so Chai-1 predicts the complex rather than the peptide alone")
    args = ap.parse_args()

    seqs = [s.strip() for s in args.input.read_text().splitlines() if s.strip()]
    receptor = None
    if args.receptor_fasta:
        lines = args.receptor_fasta.read_text().splitlines()
        receptor = "".join(l.strip() for l in lines if not l.startswith(">"))

    n = write_inputs(seqs, args.outdir, args.batches, receptor)
    print(f"wrote {n} FASTA inputs across {args.batches} batch dir(s) -> {args.outdir}")


if __name__ == "__main__":
    main()
