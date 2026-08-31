# 03 Chai-1 complex prediction

`chai_lab` v0.5.2, five complex models per sequence, MSAs from the integrated ColabFold
MMseqs2 server, no structural templates.

The receptor is the Keap1 Kelch domain, residues 326-609, in `configs/keap1_kelch.fasta`.

## Input and output

`make_chai_inputs.py` writes one FASTA per peptide, receptor as chain A and peptide as
chain B, dealt round-robin into batch directories, one per GPU. `run_chai_batch.sh` folds a
batch and writes `<name>/pred.model_idx_0-4.cif`. Finished directories are skipped.

`MSA_DIR=<dir>` serves the receptor alignment from disk instead of querying the MSA server
for every complex.
