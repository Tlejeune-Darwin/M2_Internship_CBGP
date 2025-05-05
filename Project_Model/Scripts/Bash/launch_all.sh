#!/bin/bash

# Batch number
TOTAL_BATCHES=100
SIMS_PER_BATCH=10000
MAX_JOBS=10

for i in $(seq 0 $((TOTAL_BATCHES - 1))); do
    BATCH_NAME="batch_ref_${i}"
    OFFSET=$((i * SIMS_PER_BATCH))

        # While to activate a set job number
    while true; do
        active_jobs=$(squeue -u lejeunet -h -o "%j" | grep -c "^job_run_")
        if [ "$active_jobs" -lt "$MAX_JOBS" ]; then
            break
        fi
        echo "⏳ $active_jobs Jobs running. Waiting (MAX = $MAX_JOBS)..."
        sleep 20
    done

    sbatch \
      --export=ALL,BATCH_NAME=$BATCH_NAME,NUM_SIMS=$SIMS_PER_BATCH,OFFSET=$OFFSET \
      job_run_batch.slurm
done
