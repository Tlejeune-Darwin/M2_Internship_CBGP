#!/bin/bash

# Batch number
TOTAL_BATCHES=10
SIMS_PER_BATCH=1000

for i in $(seq 0 $((TOTAL_BATCHES - 1))); do
    BATCH_NAME="batch_test_${i}"
    OFFSET=$((i * SIMS_PER_BATCH))

    sbatch \
      --export=ALL,BATCH_NAME=$BATCH_NAME,NUM_SIMS=$SIMS_PER_BATCH,OFFSET=$OFFSET \
      job_run_batch.slurm
done
