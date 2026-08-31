#!/usr/bin/env python3
"""Aggregate HADDOCK3 runs and pick one representative pose per peptide.

Walks a tree of HADDOCK3 run directories, reads the final `caprieval` table of each run,
and evaluates every reported model:

  * no score cutoff and no per-length quota -- every model is carried forward
  * models with a severe steric clash are excluded
  * models must retain contacts to the Keap1 key residues
  * where one sequence yields several poses, the electrostatic energy term (E_elec)
    picks the representative

Two CSVs are written: every model with its measurements, and one representative pose
per sequence.

Usage:
    python collect_haddock_results.py --runs haddock_runs --outdir results
    python collect_haddock_results.py --self-test
"""

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

import numpy as np

# Keap1 active residues for the docking restraints (native Keap1 numbering, chain A).
KEY_RESIDUES = [363, 380, 415, 483, 525, 530, 555, 602]

# Conventional distances: a 4.0 A heavy-atom contact and a 2.0 A severe clash. Both are
# exposed as flags; changing them changes which models are eligible, never the ranking.
CONTACT_CUTOFF = 4.0
CLASH_CUTOFF = 2.0
MIN_KEY_RESIDUES = 1

CAPRIEVAL_DIR = re.compile(r"^(\d+)_caprieval$")
RUN_SUFFIX = re.compile(r"_(\d+)$")

FIELDS = ["sequence", "length", "model_idx", "run", "model", "caprieval_rank", "score",
          "elec", "vdw", "desolv", "bsa", "dockq", "irmsd", "cluster_id",
          "min_interchain_dist", "clash", "key_residues_contacted", "n_key_residues",
          "retains_key_residues", "eligible"]


# ------------------------------------------------------------------------ structures

def read_model(path):
    """Heavy-atom coordinates of a HADDOCK model, split by chain.

    Returns (coords_A, coords_B, resids_A). Hydrogens are skipped so the distances do
    not depend on how the model was protonated.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    a_xyz, b_xyz, a_res = [], [], []

    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            element = line[76:78].strip() or line[12:16].strip()[:1]
            if element == "H":
                continue
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            chain = line[21]
            if chain == "A":
                a_xyz.append(xyz)
                a_res.append(int(line[22:26]))
            elif chain == "B":
                b_xyz.append(xyz)

    return np.array(a_xyz), np.array(b_xyz), np.array(a_res)


def measure(a_xyz, b_xyz, a_res, key_residues, contact_cutoff):
    """Minimum interchain distance and which key residues contact the peptide."""
    if len(a_xyz) == 0 or len(b_xyz) == 0:
        return float("inf"), []

    dists = np.linalg.norm(a_xyz[:, None, :] - b_xyz[None, :, :], axis=-1)
    per_atom_min = dists.min(axis=1)

    contacted = sorted({int(r) for r, d in zip(a_res, per_atom_min)
                        if d <= contact_cutoff and int(r) in set(key_residues)})
    return float(per_atom_min.min()), contacted


# ------------------------------------------------------------------------- capri I/O

def final_caprieval(run_dir):
    """Path to the highest-numbered caprieval table in a finished run, or None.

    The workflow runs caprieval several times; the last one follows seletopclusts and
    is the table the selection was made on.
    """
    candidates = []
    for path in run_dir.rglob("*_caprieval/capri_ss.tsv"):
        if "analysis" in path.parts:            # the analysis/ copies are duplicates
            continue
        m = CAPRIEVAL_DIR.match(path.parent.name)
        if m:
            candidates.append((int(m.group(1)), path))
    return max(candidates)[1] if candidates else None


def read_capri(path):
    """Rows of a capri_ss.tsv as dicts, with '-' turned into None."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [{k: (None if v == "-" else v) for k, v in row.items()} for row in rows]


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------------------------ evaluation

def evaluate_run(run_dir, key_residues, contact_cutoff, clash_cutoff, min_keys):
    """Every model of one run, measured and marked eligible or not."""
    capri = final_caprieval(run_dir)
    if capri is None:
        return []

    name = run_dir.name
    suffix = RUN_SUFFIX.search(name)
    sequence = name[:suffix.start()] if suffix else name
    model_idx = suffix.group(1) if suffix else ""

    rows = []
    for entry in read_capri(capri):
        model_path = (capri.parent / entry["model"]).resolve()
        if not model_path.exists():
            model_path = model_path.with_suffix(model_path.suffix + ".gz")
        if not model_path.exists():
            continue

        a_xyz, b_xyz, a_res = read_model(model_path)
        min_dist, contacted = measure(a_xyz, b_xyz, a_res, key_residues, contact_cutoff)

        clash = min_dist < clash_cutoff
        retains = len(contacted) >= min_keys

        rows.append({
            "sequence": sequence,
            "length": len(sequence),
            "model_idx": model_idx,
            "run": name,
            "model": Path(entry["model"]).name,
            "caprieval_rank": entry.get("caprieval_rank"),
            "score": as_float(entry.get("score")),
            "elec": as_float(entry.get("elec")),
            "vdw": as_float(entry.get("vdw")),
            "desolv": as_float(entry.get("desolv")),
            "bsa": as_float(entry.get("bsa")),
            "dockq": as_float(entry.get("dockq")),
            "irmsd": as_float(entry.get("irmsd")),
            "cluster_id": entry.get("cluster_id"),
            "min_interchain_dist": round(min_dist, 3),
            "clash": int(clash),
            "key_residues_contacted": ";".join(str(r) for r in contacted),
            "n_key_residues": len(contacted),
            "retains_key_residues": int(retains),
            "eligible": int(retains and not clash),
        })

    return rows


def select_representatives(rows):
    """One pose per sequence: the eligible model with the lowest E_elec.

    E_elec is the secondary criterion; no score cutoff and no per-length quota is
    applied, so every sequence with at least one eligible pose is represented once.
    """
    best = {}
    for row in rows:
        if not row["eligible"] or row["elec"] is None:
            continue
        current = best.get(row["sequence"])
        if current is None or row["elec"] < current["elec"]:
            best[row["sequence"]] = row
    return [best[seq] for seq in sorted(best, key=lambda s: (len(s), s))]


def self_test():
    """Check the geometry, the eligibility rules and the E_elec tie-break."""
    # chain A residue 380 sits 3 A from the peptide; residue 999 is far away
    a_xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 50.0]])
    a_res = np.array([380, 999])
    b_xyz = np.array([[0.0, 0.0, 3.0]])

    min_dist, contacted = measure(a_xyz, b_xyz, a_res, KEY_RESIDUES, 4.0)
    assert min_dist == 3.0 and contacted == [380], (min_dist, contacted)

    # a tighter contact cutoff drops the contact but not the distance
    assert measure(a_xyz, b_xyz, a_res, KEY_RESIDUES, 2.5) == (3.0, [])

    # a residue that is close but not a key residue never counts
    assert measure(a_xyz, b_xyz, np.array([381, 999]), KEY_RESIDUES, 4.0)[1] == []

    rows = [
        {"sequence": "AAA", "elec": -100.0, "eligible": 1},
        {"sequence": "AAA", "elec": -150.0, "eligible": 1},   # better E_elec -> chosen
        {"sequence": "AAA", "elec": -900.0, "eligible": 0},   # clashing -> ignored
        {"sequence": "BBBB", "elec": None, "eligible": 1},    # no energy -> ignored
        {"sequence": "CC", "elec": -10.0, "eligible": 1},
    ]
    reps = select_representatives(rows)
    assert [r["sequence"] for r in reps] == ["CC", "AAA"], reps   # sorted by length
    assert reps[1]["elec"] == -150.0, reps[1]

    assert as_float("-77.4") == -77.4 and as_float("-") is None and as_float(None) is None
    print("self-test OK")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", type=Path, help="tree containing HADDOCK3 run directories")
    ap.add_argument("--outdir", type=Path, help="where the two CSVs are written")
    ap.add_argument("--contact-cutoff", type=float, default=CONTACT_CUTOFF,
                    help="heavy-atom distance defining a key-residue contact (A)")
    ap.add_argument("--clash-cutoff", type=float, default=CLASH_CUTOFF,
                    help="interchain heavy-atom distance below which a model is a clash (A)")
    ap.add_argument("--min-key-residues", type=int, default=MIN_KEY_RESIDUES,
                    help="key-residue contacts a model must retain")
    ap.add_argument("--key-residues", type=int, nargs="+", default=KEY_RESIDUES)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return
    if not (args.runs and args.outdir):
        ap.error("--runs and --outdir are required")

    # a finished run is any directory holding numbered <n>_caprieval stage directories
    run_dirs = sorted({p.parents[1] for p in args.runs.rglob("*_caprieval/capri_ss.tsv")
                       if "analysis" not in p.parts})
    print(f"{len(run_dirs)} finished runs under {args.runs}")

    rows = []
    for i, run_dir in enumerate(run_dirs, 1):
        rows.extend(evaluate_run(run_dir, args.key_residues, args.contact_cutoff,
                                 args.clash_cutoff, args.min_key_residues))
        if i % 25 == 0 or i == len(run_dirs):
            print(f"  {i}/{len(run_dirs)} runs, {len(rows)} models", flush=True)

    args.outdir.mkdir(parents=True, exist_ok=True)
    all_csv = args.outdir / "haddock_all_models.csv"
    with open(all_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    reps = select_representatives(rows)
    summary_csv = args.outdir / "haddock_summary.csv"
    with open(summary_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(reps)

    sequences = {r["sequence"] for r in rows}
    clashing = sum(r["clash"] for r in rows)
    keeps = sum(r["retains_key_residues"] for r in rows)
    print(f"\nmodels evaluated          {len(rows)}  ({len(sequences)} sequences)")
    print(f"retaining key residues    {keeps}")
    print(f"with a steric clash       {clashing}")
    print(f"eligible                  {sum(r['eligible'] for r in rows)}")
    print(f"representative poses      {len(reps)}")
    print(f"-> {all_csv}\n-> {summary_csv}")


if __name__ == "__main__":
    sys.exit(main())
