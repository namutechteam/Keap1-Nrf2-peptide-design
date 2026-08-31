#!/bin/bash
# Run every HADDOCK3 .cfg found under a directory tree.
# Serial by default; pass a job count as the second argument to run in parallel.
#
#   ./run_haddock_batch.sh haddock_runs        # serial
#   ./run_haddock_batch.sh haddock_runs 6      # 6 concurrent runs
set -uo pipefail

ROOT_DIR="${1:?usage: $0 <root_dir> [max_parallel_jobs]}"
MAX_JOBS="${2:-1}"
LOG_DIR="$(pwd)/haddock_logs"
LOG_FILE="$LOG_DIR/haddock_run.log"
ERR_FILE="$LOG_DIR/haddock_error.log"
FAIL_FILE="$LOG_DIR/failed_jobs.txt"

mkdir -p "$LOG_DIR"
: > "$FAIL_FILE"

mapfile -t cfg_files < <(find "$ROOT_DIR" -type f -name "*.cfg" ! -path "*/analysis/*" ! -name "params.cfg" | sort)
total=${#cfg_files[@]}

echo "[$(date)] Starting HADDOCK3 batch: $total jobs, max $MAX_JOBS parallel" | tee -a "$LOG_FILE"

run_cfg() {
    local cfg="$1" idx="$2"
    local cfg_name cfg_dir
    cfg_name=$(basename "$cfg")
    cfg_dir=$(dirname "$cfg")

    echo "[$(date)] [$idx/$total] Running $cfg_name in $cfg_dir" | tee -a "$LOG_FILE"

    ( cd "$cfg_dir" && haddock3 "$cfg_name" ) >> "$LOG_FILE" 2>> "$ERR_FILE"

    if [ $? -eq 0 ]; then
        echo "[OK]   $cfg_name" | tee -a "$LOG_FILE"
    else
        echo "[FAIL] $cfg_name" | tee -a "$LOG_FILE" "$ERR_FILE"
        echo "$cfg" >> "$FAIL_FILE"
    fi
}

count=0
for cfg in "${cfg_files[@]}"; do
    count=$((count + 1))
    if [ "$MAX_JOBS" -le 1 ]; then
        run_cfg "$cfg" "$count"
    else
        run_cfg "$cfg" "$count" &
        while [ "$(jobs -r | wc -l)" -ge "$MAX_JOBS" ]; do sleep 1; done
    fi
done
wait

fail_count=$(wc -l < "$FAIL_FILE")
echo "[$(date)] Finished. Success: $((total - fail_count))/$total | Failed: $fail_count" | tee -a "$LOG_FILE"
