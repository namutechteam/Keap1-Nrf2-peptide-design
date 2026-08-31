# 05 HADDOCK3 refinement

HADDOCK3 v2024.10.0b7.

| | |
|---|---|
| receptor | Keap1 Kelch domain, chain A, residues 326-609 |
| restraints | seven ambiguous interaction restraints, `2.0 2.0 0.0` |
| Keap1 active residues | Ser363, Arg380, Arg415, Arg483, Tyr525, Gln530, Ser555, Ser602 |
| peptide active residues | Nrf2 78, 79, 80 and 82 |
| rigid-body models | 2,000 per complex |
| carried to flexref | top 100 |
| refinement | flexref then emref, peptide chain B fully flexible |
| clustering | fraction of common contacts |

`make_air_tbl.py` writes the restraint file for any peptide length from one definition in
Nrf2 residue numbering.

## Scripts

`make_haddock_runs.py` builds one run directory per peptide PDB, reading the length from
the file. `run_haddock_batch.sh` runs the jobs. `collect_haddock_results.py` reads the
final `caprieval` table of each run, measures clashes and contacts to the restrained
residues, and picks one representative pose per sequence by electrostatic energy.
