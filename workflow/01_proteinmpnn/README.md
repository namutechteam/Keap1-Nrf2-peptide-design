# 01 ProteinMPNN sequence generation

ProteinMPNN v1.0.1, vanilla full-backbone model `v_48_020`.

| | |
|---|---|
| template | PDB 4IFL, Keap1 chain X held fixed, peptide chain P designed |
| peptide lengths | 10-16 residues, by sequential N-terminal truncation |
| sampling temperature | 0.5 |
| random seed | 42 |
| sequences per length | 100,000 |
| batch_size | 16 for the 10-15mers, 100 for the 16mer |

## Scripts

`prepare_inputs.sh` builds the two JSONL files from a trimmed template.
`run_proteinmpnn.sh` generates the designs for one length. Set `MPNN_DIR` to a ProteinMPNN
checkout; `BATCH_SIZE` overrides the default.

## Output

ProteinMPNN FASTA, `seqs/4IFL_<L>.fa`. The first record is the native template sequence,
each following record is one design.
