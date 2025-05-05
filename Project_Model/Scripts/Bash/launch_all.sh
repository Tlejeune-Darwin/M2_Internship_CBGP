#!/bin/bash

# Nombre de batchs
TOTAL_BATCHES=5
SIMS_PER_BATCH=10
MAX_JOBS=5

for i in $(seq 0 $((TOTAL_BATCHES - 1))); do
    BATCH_NAME="batch_${i}"
    OFFSET=$((i * SIMS_PER_BATCH))

    echo "Lancement de $BATCH_NAME avec OFFSET=$OFFSET"

    sbatch \
      --export=ALL,BATCH_NAME=$BATCH_NAME,NUM_SIMS=$SIMS_PER_BATCH,OFFSET=$OFFSET \
      job_run_batch.slurm
done
