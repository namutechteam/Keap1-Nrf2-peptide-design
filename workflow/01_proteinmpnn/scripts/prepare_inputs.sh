#!/bin/bash
# Build the ProteinMPNN input JSONL for one peptide length from the trimmed 4IFL template.
#
# Uses ProteinMPNN's own helper scripts, so point MPNN_DIR at your ProteinMPNN checkout
# (commit 8907e6671bfbfc92303b5f79c4b5e6ce47cdef57).
#
#   ./prepare_inputs.sh 12 ../inputs/4IFL_12.pdb ../inputs
#   -> ../inputs/12mer.jsonl and ../inputs/12mer_P_chain.jsonl
set -euo pipefail

LENGTH="${1:?usage: $0 <peptide_length> <template_pdb> <outdir>}"
TEMPLATE="${2:?usage: $0 <peptide_length> <template_pdb> <outdir>}"
OUT_DIR="${3:?usage: $0 <peptide_length> <template_pdb> <outdir>}"
MPNN_DIR="${MPNN_DIR:-$HOME/ProteinMPNN}"

mkdir -p "$OUT_DIR"
PDB_DIR=$(mktemp -d)
trap 'rm -rf "$PDB_DIR"' EXIT
cp "$TEMPLATE" "$PDB_DIR/"

# chain X = Keap1 Kelch domain (fixed), chain P = peptide (designed)
python "$MPNN_DIR/helper_scripts/parse_multiple_chains.py" \
    --input_path "$PDB_DIR" \
    --output_path "$OUT_DIR/${LENGTH}mer.jsonl"

python "$MPNN_DIR/helper_scripts/assign_fixed_chains.py" \
    --input_path "$OUT_DIR/${LENGTH}mer.jsonl" \
    --output_path "$OUT_DIR/${LENGTH}mer_P_chain.jsonl" \
    --chain_list "P"

echo "wrote $OUT_DIR/${LENGTH}mer.jsonl and $OUT_DIR/${LENGTH}mer_P_chain.jsonl"
