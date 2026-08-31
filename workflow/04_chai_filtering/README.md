# 04 Structure-based filtering

Two criteria, both required, defined in `configs/structure_filters.json`.

| | |
|---|---|
| intrapeptide hydrogen bond | at least one backbone N-O contact within 3.0 A, ignoring the donor residue and its two neighbours |
| Keap1 hot spots | a hydrogen bond to Ser363, Arg380, Asn382, Arg415, Arg483 and Ser508 |

PLIP runs in peptide-chain interaction mode. Hot spots are listed in Keap1 numbering and
converted with `receptor_start`.

## Input and output

`structure_filter.py` reads a tree of predicted `.cif` files at any depth and writes one
row per model to `processed_summary.csv`: `sequence, model_idx, cif_file, pdb_file,
intrahbond_count, ppi_profile, hotspots_missing, pass, error`.

`cif_to_haddock_peptide.py` writes the peptide of each passing model as
`<SEQUENCE>_<model_idx>.pdb` for stage 05.
