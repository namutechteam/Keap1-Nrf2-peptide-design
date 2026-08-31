# analysis

Post-hoc analysis. These scripts read `data/` and write the manuscript tables and
figures.

Run in the `keap1nrf2` environment created by `environment/setup.sh`.

| script | produces |
|---|---|
| `cascade_composition_audit.py` | Table 4, Table S4, Table S5 |
| `haddock_checkpoint_composition.py` | the 10-15-residue block of Table S4 |
| `retrospective_metrics.py` | every computed quantity against measured inhibition |
| `build_dataset.py` | rebuilds `Keap1_Nrf2_45peptide_dataset.csv` from `screening_table.csv` and the Chai-1 scores |
| `make_figures.py` | Figure 4, and the computed quantities against measured inhibition |
| `make_toc.py` | the table-of-contents graphic |

## Running them

```bash
python analysis/code/cascade_composition_audit.py \
       analysis/data/sequences analysis/data/Keap1_Nrf2_45peptide_dataset.csv analysis/data/derived

python analysis/code/haddock_checkpoint_composition.py \
       analysis/data/sequences analysis/data/sequences/03_haddock3_checkpoint_10to15mer.csv \
       analysis/data/Keap1_Nrf2_45peptide_dataset.csv \
       analysis/data/derived/si_table_S4_haddock.json

python analysis/code/retrospective_metrics.py \
       analysis/data/Keap1_Nrf2_45peptide_dataset.csv analysis/data/derived/table_metric_activity.csv

python analysis/code/build_dataset.py \
       analysis/data/screening_table.csv \
       analysis/data/chai1_confidence/scores_all_models.csv \
       rebuilt_dataset.csv

# make_figures.py reads si_table_S4_haddock.json, so run the step above first
python analysis/code/make_figures.py analysis/data/derived analysis/data/Keap1_Nrf2_45peptide_dataset.csv figures
python analysis/code/make_toc.py figures
```

`cascade_composition_audit.py` and `haddock_checkpoint_composition.py` take the
sequence directory as their first argument and read the `01_generated/<L>mer/` and
`02_filtered/<L>mer/` layout used in `data/sequences`.
`haddock_checkpoint_composition.py` reads the docking checkpoint either from the CSV
above or from the original workbook, whichever is passed.

## What these scripts do not cover

The selection of 28 candidates from the 140 docking models included a visual
review step. Its criteria are stated in the Methods and in Supporting Information
1, but they were not reduced to a numerical threshold, so that step is not
reproduced here and the final selection is not exactly recoverable from these
files. This is stated in the paper.
