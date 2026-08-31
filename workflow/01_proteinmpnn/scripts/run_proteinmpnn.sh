#!/bin/bash
# Generate the design set for one peptide length.
#
#   ./run_proteinmpnn.sh 12 ../inputs ../outputs
#
# Point MPNN_DIR at your ProteinMPNN checkout.
set -euo pipefail

LENGTH="${1:?usage: $0 <peptide_length> <input_dir> <output_dir>}"
IN_DIR="${2:?usage: $0 <peptide_length> <input_dir> <output_dir>}"
OUT_DIR="${3:?usage: $0 <peptide_length> <input_dir> <output_dir>}"
MPNN_DIR="${MPNN_DIR:-$HOME/ProteinMPNN}"

# batch_size is not a free performance knob: protein_mpnn_run.py draws
# torch.randn(batch_size, L) once per batch, so it determines the whole sampling
# stream. The same seed with a different batch_size gives a different sequence set.
if [ "$LENGTH" -ge 16 ]; then BATCH_SIZE="${BATCH_SIZE:-100}"; else BATCH_SIZE="${BATCH_SIZE:-16}"; fi

mkdir -p "$OUT_DIR"

python "$MPNN_DIR/protein_mpnn_run.py" \
    --jsonl_path        "$IN_DIR/${LENGTH}mer.jsonl" \
    --chain_id_jsonl    "$IN_DIR/${LENGTH}mer_P_chain.jsonl" \
    --out_folder        "$OUT_DIR" \
    --model_name        v_48_020 \
    --num_seq_per_target 100000 \
    --sampling_temp     0.5 \
    --batch_size        "$BATCH_SIZE" \
    --seed              42 \
    2>&1 | tee "$OUT_DIR/run_${LENGTH}mer.log"

echo "designs -> $OUT_DIR/seqs/4IFL_${LENGTH}.fa"
