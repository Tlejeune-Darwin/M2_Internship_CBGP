#!/bin/bash

# Nombre de batchs
TOTAL_BATCHES=4
SIMS_PER_BATCH=10
MAX_JOBS=1

for i in $(seq 0 $((TOTAL_BATCHES - 1))); do
    BATCH_NAME="batch_ref_${i}"
    OFFSET=$((i * SIMS_PER_BATCH))

        # Tant qu'on a déjà trop de jobs en attente ou en cours...
    while [ $(squeue -u lejeunet -h -o "%j" | grep -c "^batch_") -ge $MAX_JOBS ]; do
        echo "⏳ $MAX_JOBS jobs déjà en cours. Attente..."
        sleep 20
    done

    echo "Lancement de $BATCH_NAME avec OFFSET=$OFFSET"

    sbatch \
      --export=ALL,BATCH_NAME=$BATCH_NAME,NUM_SIMS=$SIMS_PER_BATCH,OFFSET=$OFFSET \
      job_run_batch.slurm
done
