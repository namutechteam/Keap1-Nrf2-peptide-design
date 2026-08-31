#!/bin/bash
# Run Chai-1 over one batch directory of FASTA inputs.
# Launch once per GPU:
#   ./run_chai_batch.sh chai_in/batch0 chai_out cuda:0 &
#   ./run_chai_batch.sh chai_in/batch1 chai_out cuda:1 &
set -euo pipefail

IN_DIR="${1:?usage: $0 <input_dir> <output_dir> [device]}"
OUT_DIR="${2:?usage: $0 <input_dir> <output_dir> [device]}"
DEVICE="${3:-cuda:0}"

# Optional: serve the receptor alignment from disk instead of querying ColabFold for
# every complex. The receptor is the same query each time, so the MSA context is
# identical either way -- see workflow/03_chai1/README.md. Peptides return no hits and fall back
# to single-sequence, which is what the server returns for them anyway.
#   MSA_DIR=shared_msa ./run_chai_batch.sh chai_in/batch0 chai_out cuda:0
MSA_DIR="${MSA_DIR:-}"
if [ -n "$MSA_DIR" ]; then
    MSA_ARGS=(--msa-directory "$MSA_DIR")
else
    MSA_ARGS=(--use-msa-server)
fi

mkdir -p "$OUT_DIR"

for fasta in "$IN_DIR"/*.fasta; do
    name=$(basename "$fasta" .fasta)
    dest="$OUT_DIR/$name"

    # skip work already done, so an interrupted run can just be restarted
    if [ -d "$dest" ] && compgen -G "$dest/*.cif" > /dev/null; then
        continue
    fi

    # chai fold asserts its output directory is empty; an interrupted run leaves one
    # holding only msas/, which would fail every retry, so clear it first
    rm -rf "$dest"
    mkdir -p "$dest"
    chai fold "${MSA_ARGS[@]}" --device "$DEVICE" "$fasta" "$dest"
done

echo "[$(date)] finished $IN_DIR on $DEVICE"
