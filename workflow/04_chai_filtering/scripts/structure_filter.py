#!/usr/bin/env python3
"""Structure-based filtering of Chai-1 predictions.

For every predicted model this computes two independent criteria and writes one row
per model to `processed_summary.csv`:

1. intra-H-bond filter: the peptide backbone must form at least one N...O contact within
   intra_hbond_cutoff (3.0 A), ignoring the donor's own residue and its two sequence
   neighbours so that peptide bonds do not count.

2. hotspot filter: PLIP is run in peptide mode and the model is kept only if it forms a
   hydrogen bond to every Keap1 hotspot residue listed in the config.

A model passes step 4 when both criteria hold.

Hotspots are written in native Keap1 residue numbers in the config and converted to the
prediction's own numbering with `receptor_start`, so the config stays readable and does
not have to be rewritten if the receptor construct changes.

Usage:
    python structure_filter.py --input chai_out --outdir filter_out \
                               --config ../inputs/structure_filters.json --jobs 8
    python structure_filter.py --self-test
"""

import argparse
import csv
import json
import math
import re
import sys
from functools import partial
from multiprocessing import Pool
from pathlib import Path

THREE2ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
    "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
    "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
    "TRP": "W", "TYR": "Y", "ASX": "B", "GLX": "Z", "XLE": "J", "XAA": "X",
}

MODEL_IDX = re.compile(r"model_idx_(\d+)")


# --------------------------------------------------------------------------- PDB I/O

def parse_backbone(pdb_text, chain):
    """Backbone N and O atoms of one chain as (serial, atom_name, resseq, x, y, z).

    Parsed by PDB column positions rather than by splitting on whitespace: coordinates
    of large structures run together and would otherwise be mis-assigned.
    """
    atoms = []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM") or line[21] != chain:
            continue
        name = line[12:16].strip()
        if name not in ("N", "O"):              # backbone amide N / carbonyl O only
            continue
        atoms.append((int(line[6:11]), name, int(line[22:26]),
                      float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return atoms


def intra_hbond_pairs(atoms, cutoff):
    """N...O pairs within cutoff, skipping the donor's own residue and its neighbours."""
    donors = [a for a in atoms if a[1] == "N"]
    acceptors = [a for a in atoms if a[1] == "O"]

    pairs = []
    for serial, _, res, x, y, z in donors:
        for o_serial, _, o_res, ox, oy, oz in acceptors:
            if o_res in (res - 1, res, res + 1):
                continue
            d = math.dist((x, y, z), (ox, oy, oz))
            if d <= cutoff:
                pairs.append((serial, o_serial, round(d, 3)))
    return pairs


# ------------------------------------------------------------------------- one model

def profile_model(cif_path, outdir, chain, cutoff, hotspots):
    """Convert, measure and profile a single Chai-1 CIF. Returns one summary row."""
    from Bio.PDB import PDBIO
    from Bio.PDB.MMCIFParser import MMCIFParser
    from plip.basic import config as plip_config
    from plip.exchange.report import BindingSiteReport
    from plip.structure.preparation import PDBComplex

    cif_path = Path(cif_path)
    idx = MODEL_IDX.search(cif_path.name)
    model_idx = int(idx.group(1)) if idx else -1

    structure = MMCIFParser(QUIET=True).get_structure("m", str(cif_path))
    peptide = structure[0][chain]
    sequence = "".join(THREE2ONE.get(r.resname, "X") for r in peptide.get_residues())

    pdb_path = Path(outdir) / f"{sequence}_{model_idx}.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdb_path))

    pairs = intra_hbond_pairs(parse_backbone(pdb_path.read_text(), chain), cutoff)

    # PLIP in peptide mode: the peptide chain is characterised as the ligand.
    plip_config.PEPTIDES = [chain]
    complex_ = PDBComplex()
    complex_.load_pdb(str(pdb_path))
    complex_.analyze()
    for ligand in complex_.ligands:
        complex_.characterize_complex(ligand)

    profile = set()
    for _, site in sorted(complex_.interaction_sets.items()):
        report = BindingSiteReport(site)
        for kind in ("hbond", "hydrophobic", "saltbridge"):
            columns = getattr(report, f"{kind}_features")
            resnr = columns.index("RESNR")
            for row in getattr(report, f"{kind}_info"):
                profile.add(f"{row[resnr]}_{kind}")
        break                                    # the peptide gives a single binding site

    return {
        "sequence": sequence,
        "model_idx": model_idx,
        "cif_file": cif_path.name,
        "pdb_file": pdb_path.name,
        "intrahbond_count": len(pairs),
        "ppi_profile": ";".join(sorted(profile)),
        "hotspots_missing": ";".join(sorted(hotspots - profile)),
        "pass": int(bool(pairs) and hotspots <= profile),
    }


def _safe_profile(cif_path, **kwargs):
    try:
        return profile_model(cif_path, **kwargs)
    except Exception as exc:                     # one bad structure must not stop the sweep
        return {"sequence": "", "model_idx": -1, "cif_file": Path(cif_path).name,
                "pdb_file": "", "intrahbond_count": 0, "ppi_profile": "",
                "hotspots_missing": "", "pass": 0, "error": f"{type(exc).__name__}: {exc}"}


# ---------------------------------------------------------------------------- config

def hotspot_labels(config):
    """Config hotspots (native Keap1 numbering) -> PLIP labels in prediction numbering.

    Chai-1 renumbers the receptor from 1, so a hotspot at Keap1 residue R appears as
    R - receptor_start + 1.
    """
    offset = config["receptor_start"] - 1
    return {f"{h['residue'] - offset}_{h['interaction']}" for h in config["hotspots"]}


FIELDS = ["sequence", "model_idx", "cif_file", "pdb_file", "intrahbond_count",
          "ppi_profile", "hotspots_missing", "pass", "error"]


def self_test():
    """Check the geometry, numbering and pass/fail logic without PLIP or Biopython."""
    pdb = "\n".join([
        # chain B backbone: N of residue 5 sits 2.5 A from O of residue 1 -> one H-bond
        "ATOM      1  N   GLY B   1       0.000   0.000   0.000  1.00  0.00           N",
        "ATOM      2  O   GLY B   1       0.000   0.000  10.000  1.00  0.00           O",
        "ATOM      3  N   GLY B   2       0.000   0.000   3.000  1.00  0.00           N",
        "ATOM      4  O   GLY B   2       0.000   0.000   4.000  1.00  0.00           O",
        "ATOM      5  N   GLY B   5       0.000   0.000  12.500  1.00  0.00           N",
        "ATOM      6  CA  GLY B   5       0.000   0.000  13.000  1.00  0.00           C",
        "ATOM      7  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N",
    ])
    atoms = parse_backbone(pdb, "B")
    assert len(atoms) == 5, atoms                     # chain A and the CA are excluded

    pairs = intra_hbond_pairs(atoms, 3.0)
    # N(res5)-O(res1) = 2.5 A qualifies; N(res2)-O(res1) and N(res2)-O(res2) are neighbours
    assert pairs == [(5, 2, 2.5)], pairs
    assert intra_hbond_pairs(atoms, 2.0) == []        # cutoff is honoured

    # Keap1 Ser363 with receptor construct starting at 326 -> PLIP index 38
    cfg = {"receptor_start": 326,
           "hotspots": [{"residue": 363, "interaction": "hbond"},
                        {"residue": 508, "interaction": "hbond"}]}
    labels = hotspot_labels(cfg)
    assert labels == {"38_hbond", "183_hbond"}, labels

    # pass requires an intra H-bond AND every hotspot present
    assert labels <= {"38_hbond", "183_hbond", "90_hbond"}
    assert not labels <= {"38_hbond"}
    print("self-test OK")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, help="directory of per-peptide Chai-1 outputs")
    ap.add_argument("--outdir", type=Path, help="converted PDBs and processed_summary.csv")
    ap.add_argument("--config", type=Path, help="JSON hotspot / cutoff definition")
    ap.add_argument("--jobs", type=int, default=1, help="worker processes")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return
    if not (args.input and args.outdir and args.config):
        ap.error("--input, --outdir and --config are required")

    config = json.loads(args.config.read_text())
    hotspots = hotspot_labels(config)
    args.outdir.mkdir(parents=True, exist_ok=True)

    # search recursively, so predictions may sit directly under --input or be
    # grouped by peptide length
    cifs = sorted(str(p) for p in args.input.rglob("*.cif"))
    print(f"{len(cifs)} models; requiring hotspots {sorted(hotspots)}")

    work = partial(_safe_profile, outdir=str(args.outdir), chain=config["peptide_chain"],
                   cutoff=config["intra_hbond_cutoff"], hotspots=hotspots)
    if args.jobs > 1:
        with Pool(args.jobs) as pool:
            rows = pool.map(work, cifs)
    else:
        rows = [work(c) for c in cifs]

    summary = args.outdir / "processed_summary.csv"
    with open(summary, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    errors = sum("error" in r for r in rows)
    with_hbond = sum(r["intrahbond_count"] > 0 for r in rows)
    passed = sum(r["pass"] for r in rows)
    print(f"\nmodels processed        {len(rows)}  ({errors} failed)")
    print(f"after intra-H-bond      {with_hbond}")
    print(f"after hotspot filter    {passed}")
    print(f"-> {summary}")


if __name__ == "__main__":
    sys.exit(main())
