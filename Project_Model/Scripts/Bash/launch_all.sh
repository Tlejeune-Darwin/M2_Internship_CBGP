#!/bin/bash

# Nombre de batchs
TOTAL_BATCHES=100
SIMS_PER_BATCH=10000
MAX_JOBS=10

for i in $(seq 0 $((TOTAL_BATCHES - 1))); do
    BATCH_NAME="batch_${i}"
    OFFSET=$((i * SIMS_PER_BATCH))

        # Tant qu'on a déjà trop de jobs en attente ou en cours...
    while [ $(squeue -u $USER -h -n $BATCH_NAME | wc -l) -ge $MAX_JOBS ]; do
        echo "⏳ $MAX_JOBS jobs déjà en cours. Attente..."
        sleep 20
    done

    echo "Lancement de $BATCH_NAME avec OFFSET=$OFFSET"

    sbatch \
      --export=ALL,BATCH_NAME=$BATCH_NAME,NUM_SIMS=$SIMS_PER_BATCH,OFFSET=$OFFSET \
      job_run_batch.slurm
done
