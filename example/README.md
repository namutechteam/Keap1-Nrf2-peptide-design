# Example

The 12-residue design run. Every stage has its input here, so any stage can be run on its
own without redoing the one before it.

```
01_proteinmpnn/4IFL_12.pdb            template trimmed to 12 residues
01_proteinmpnn/12mer.jsonl            parsed template
01_proteinmpnn/12mer_P_chain.jsonl    chain assignment, P designed and X fixed
02_sequence_filtering/4IFL_12.fa      100,000 generated sequences
03_chai1/filtered_12mer.txt           the 770 sequences that passed the filters
04_chai_filtering/chai_out/           predicted complexes for two of them
05_haddock3/12mer.tbl                 AIR restraints for a 12mer
05_haddock3/protein.pdb               prepared Keap1 receptor, chain A
05_haddock3/reference.pdb             4IFL complex, caprieval reference
```

Run from the repository root, with the environments from `environment/setup.sh`.

```bash
conda activate keap1nrf2

# 01  generate sequences, about 3 h on one GPU
export MPNN_DIR=~/ProteinMPNN
workflow/01_proteinmpnn/scripts/run_proteinmpnn.sh 12 example/01_proteinmpnn out_12mer

# 02  sequence filtering, seconds
python workflow/02_sequence_filtering/scripts/filter_sequences.py \
    --input   example/02_sequence_filtering/4IFL_12.fa \
    --config  workflow/02_sequence_filtering/configs/keap1_filters.json \
    --output  filtered_12.txt \
    --summary counts_12.csv

# 03  fold the survivors, about a day on two GPUs
python workflow/03_chai1/scripts/make_chai_inputs.py \
    --input example/03_chai1/filtered_12mer.txt --outdir chai_in --batches 2 \
    --receptor-fasta workflow/03_chai1/configs/keap1_kelch.fasta
workflow/03_chai1/scripts/run_chai_batch.sh chai_in/batch0 chai_out cuda:0

# 04  structure filtering, seconds on the two supplied predictions
python workflow/04_chai_filtering/scripts/structure_filter.py \
    --input example/04_chai_filtering/chai_out --outdir filter_out \
    --config workflow/04_chai_filtering/configs/structure_filters.json --jobs 8
python workflow/04_chai_filtering/scripts/cif_to_haddock_peptide.py \
    --input example/04_chai_filtering/chai_out --outdir haddock_input \
    --passed filter_out/processed_summary.csv

# 05  docking, about 1.7 h per run
conda activate keap1nrf2-haddock
python workflow/05_haddock3/scripts/make_haddock_runs.py \
    --peptides  haddock_input \
    --receptor  example/05_haddock3/protein.pdb \
    --reference example/05_haddock3/reference.pdb \
    --template  workflow/05_haddock3/configs/haddock3_template.cfg \
    --outdir    runs
workflow/05_haddock3/scripts/run_haddock_batch.sh runs 6
python workflow/05_haddock3/scripts/collect_haddock_results.py --runs runs --outdir results
```

Stages 02, 04 and 05 run from the files here. Stages 01 and 03 are the expensive ones;
their outputs are supplied, so they can be skipped.

`make_air_tbl.py` writes the restraint file for any other peptide length.
