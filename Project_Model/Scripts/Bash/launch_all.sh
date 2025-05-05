#!/bin/bash

# Nombre de batchs
TOTAL_BATCHES=2
SIMS_PER_BATCH=10
MAX_JOBS=10

RESULTS_DIR="$HOME/results"

start_index=0
for i in $(seq 0 $((TOTAL_BATCHES - 1))); do
    if [ ! -d "$RESULTS_DIR/batch_$i" ]; then
        start_index=$i
        break
    fi
done
echo "⏩ Reprise à partir de batch_$start_index"

for i in $(seq $start_index $((TOTAL_BATCHES - 1))); do
    BATCH_NAME="batch_${i}"
    OFFSET=$((i * SIMS_PER_BATCH))

    echo "🚀 Lancement de $BATCH_NAME avec OFFSET=$OFFSET"

    sbatch \
      --export=ALL,BATCH_NAME=$BATCH_NAME,NUM_SIMS=$SIMS_PER_BATCH,OFFSET=$OFFSET \
      job_run_batch.slurm
done
