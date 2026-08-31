#!/usr/bin/env python3
"""Build one self-contained HADDOCK3 run directory per peptide ensemble.

Each run directory gets the receptor, the reference complex, the peptide
ensemble, a length-matched AIR restraint file (from make_air_tbl.py) and a .cfg
rendered from the template with the peptide name and length substituted in.

Usage:
    python make_haddock_runs.py \
        --peptides chai_out/result --receptor protein.pdb --reference reference.pdb \
        --template ../inputs/haddock3_template.cfg --outdir haddock_runs
"""

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from make_air_tbl import build_tbl                       # noqa: E402


def peptide_length(pdb_path):
    """Residue count of the first model of a peptide ensemble PDB."""
    resids = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM"):
            key = (line[21], line[22:27])                # chain + resSeq/icode
            if key not in resids:
                resids.append(key)
    return len(resids)


def make_run(pdb, receptor, reference, template, outdir):
    name = pdb.stem
    length = peptide_length(pdb)
    run_dir = outdir / name
    run_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy(receptor, run_dir / "protein.pdb")
    shutil.copy(reference, run_dir / "reference.pdb")
    shutil.copy(pdb, run_dir / f"{name}.pdb")
    (run_dir / "tbl.tbl").write_text(build_tbl(length))
    (run_dir / f"{name}.cfg").write_text(
        template.replace("{{PEPTIDE}}", name).replace("{{PEPTIDE_LENGTH}}", str(length)))

    return name, length


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--peptides", type=Path, required=True,
                    help="directory of peptide ensemble PDBs (output of step 04)")
    ap.add_argument("--receptor", type=Path, required=True, help="prepared Keap1 PDB")
    ap.add_argument("--reference", type=Path, required=True, help="reference complex for caprieval")
    ap.add_argument("--template", type=Path, required=True, help="haddock3 .cfg template")
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    template = args.template.read_text()
    for pdb in sorted(args.peptides.glob("*.pdb")):
        name, length = make_run(pdb, args.receptor, args.reference, template, args.outdir)
        print(f"{name}\t{length} aa")


if __name__ == "__main__":
    main()
