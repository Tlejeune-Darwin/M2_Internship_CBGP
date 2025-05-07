#!/bin/bash

# Demander à l'utilisateur les infos essentielles
read -p "Simulations file name : " SIM_DIR
read -p "Start from an existant simulation ? (y/n) : " RESUME
read -p "Batch number : " TOTAL_BATCHES
read -p "Simulation per batch : " SIMS_PER_BATCH

# Chemin de base vers le dossier de résultats
BASE_PATH="$HOME/results/$SIM_DIR"
mkdir -p "$BASE_PATH"

# Déterminer le point de départ
if [[ "$RESUME" == "y" ]]; then
    # Compte les batchs déjà présents
    LAST_BATCH=$(ls "$BASE_PATH" | grep "^batch_${SIM_DIR}_" | wc -l)
    START_BATCH=$LAST_BATCH
    echo "Starting from batch $START_BATCH"
else
    START_BATCH=0
    echo "Starting from scratch"
fi

# Lancer les batchs
for i in $(seq $START_BATCH $((START_BATCH + TOTAL_BATCHES - 1))); do
    BATCH_NAME="batch_${SIM_DIR}_${i}"
    OFFSET=$((i * SIMS_PER_BATCH))

    echo "Launching $BATCH_NAME with OFFSET=$OFFSET"

    sbatch \
      --export=ALL,BATCH_NAME=$BATCH_NAME,NUM_SIMS=$SIMS_PER_BATCH,OFFSET=$OFFSET \
      job_run_batch.slurm
done
