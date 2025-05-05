#!/bin/bash

# Nombre de batchs
TOTAL_BATCHES=4
SIMS_PER_BATCH=10000
MAX_JOBS=4

for i in $(seq 0 $((TOTAL_BATCHES - 1))); do
    BATCH_NAME="batch_${i}"
    OFFSET=$((i * SIMS_PER_BATCH))

    sbatch \
      --export=ALL,BATCH_NAME=$BATCH_NAME,NUM_SIMS=$SIMS_PER_BATCH,OFFSET=$OFFSET \
      --job-name=$BATCH_NAME \
      job_run_batch.slurm
done
