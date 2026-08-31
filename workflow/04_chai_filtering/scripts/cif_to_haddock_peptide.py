#!/usr/bin/env python3
"""Extract the peptide chain from Chai-1 CIF predictions as HADDOCK3 input PDBs.

Chai-1 predicts the Keap1-peptide complex; HADDOCK3 re-docks the peptide against a
separately prepared receptor. Each Chai model that passed the step-04 structural
filters is written out as its own single-model PDB:

    <SEQUENCE>_<model_idx>.pdb

One HADDOCK3 run is then set up per file, which is why a sequence can appear with
one or two model indices rather than all five.

Input layout (as produced by workflow/03_chai1/scripts/run_chai_batch.sh):
    chai_out/0/pred.model_idx_0.cif ... pred.model_idx_4.cif
    chai_out/1/...

Usage:
    # convert every model of every prediction
    python cif_to_haddock_peptide.py --input chai_out --outdir haddock_input

    # convert only the models listed by the structural filter (CSV: sequence,model_idx)
    python cif_to_haddock_peptide.py --input chai_out --outdir haddock_input \
                                     --passed processed_summary.csv

    python cif_to_haddock_peptide.py --self-test
"""

import argparse
import csv
import re
import sys
from pathlib import Path

THREE2ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
    "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
    "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
    "TRP": "W", "TYR": "Y", "ASX": "B", "GLX": "Z", "XLE": "J", "XAA": "X",
}

# Chai-1 names its outputs pred.model_idx_<N>.cif
MODEL_IDX = re.compile(r"model_idx_(\d+)")


def model_index(cif_path):
    m = MODEL_IDX.search(cif_path.name)
    return int(m.group(1)) if m else None


def sequence_of(chain):
    return "".join(THREE2ONE.get(r.resname, "X") for r in chain.get_residues())


def shortest_chain(model):
    """The peptide is the shorter of the two chains in a Keap1-peptide complex."""
    return min(model, key=lambda c: len(list(c.get_residues())))


def read_passed(csv_path):
    """Read the structural filter's verdict as {(sequence, model_idx)}.

    Accepts any CSV with 'sequence' and 'model_idx' columns; if a 'pass' column is
    present, only truthy rows are kept.
    """
    passed = set()
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            verdict = row.get("pass", "1").strip().lower()
            if verdict in ("0", "false", "no", "fail", ""):
                continue
            passed.add((row["sequence"].strip(), int(row["model_idx"])))
    return passed


def convert(pred_dir, outdir, chain_id=None, passed=None):
    """Convert one prediction directory. Returns the list of written paths."""
    from Bio.PDB import PDBIO, Select
    from Bio.PDB.MMCIFParser import MMCIFParser

    class ChainOnly(Select):
        def __init__(self, keep):
            self.keep = keep

        def accept_chain(self, chain):
            return chain.id == self.keep

    parser, io = MMCIFParser(QUIET=True), PDBIO()
    written = []

    for cif in sorted(pred_dir.glob("*.cif")):
        idx = model_index(cif)
        if idx is None:
            continue

        model = parser.get_structure("m", cif)[0]
        chain = model[chain_id] if chain_id else shortest_chain(model)
        seq = sequence_of(chain)

        if passed is not None and (seq, idx) not in passed:
            continue

        outdir.mkdir(parents=True, exist_ok=True)
        out = outdir / f"{seq}_{idx}.pdb"
        io.set_structure(model)
        io.save(str(out), ChainOnly(chain.id))
        written.append(out)

    return written


def self_test():
    """Check the pieces that do not need Biopython or real CIF files."""
    assert model_index(Path("pred.model_idx_3.cif")) == 3
    assert model_index(Path("pred.model_idx_0.cif")) == 0
    assert model_index(Path("something_else.cif")) is None

    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, newline="") as f:
        f.write("sequence,model_idx,pass\n"
                "TLVLDPRTGELS,1,1\n"
                "TLVLDPRTGELS,2,0\n"          # explicit fail -> excluded
                "ALTLDPHTGELL,3,true\n")
        path = f.name
    passed = read_passed(path)
    Path(path).unlink()

    assert passed == {("TLVLDPRTGELS", 1), ("ALTLDPHTGELL", 3)}, passed
    print("self-test OK")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, help="directory of per-peptide Chai-1 outputs")
    ap.add_argument("--outdir", type=Path, help="where the peptide PDBs are written")
    ap.add_argument("--chain", help="peptide chain ID (default: the shorter chain)")
    ap.add_argument("--passed", type=Path,
                    help="CSV of models that passed the structural filter "
                         "(columns: sequence, model_idx, optional pass)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return
    if not (args.input and args.outdir):
        ap.error("--input and --outdir are required")

    passed = read_passed(args.passed) if args.passed else None
    total = 0
    for pred_dir in sorted(p for p in args.input.iterdir() if p.is_dir()):
        for out in convert(pred_dir, args.outdir, args.chain, passed):
            print(out.name)
            total += 1
    print(f"\nwrote {total} peptide PDBs -> {args.outdir}")


if __name__ == "__main__":
    sys.exit(main())
