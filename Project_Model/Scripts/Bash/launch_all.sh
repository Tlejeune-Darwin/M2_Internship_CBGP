#!/bin/bash

#  Asking the user some inputs
read -p "Simulations file name : " SIM_DIR                          # The name used to identify the file in the results folder
read -p "Start from an existant simulation ? (y/n) : " RESUME       # You need to inform the script if there are simulations in the folder already. If yes, it will start from the highest simulation number
read -p "Batch number : " TOTAL_BATCHES                             # Number of batches to run, 1 batch is equal to 1 job
read -p "Simulation per batch : " SIMS_PER_BATCH                    # Number of simulation per batch, the more the simulations, the more time the batch will be running in the cluster

# Base file and creation of the folder name input by the user
BASE_PATH="$HOME/results/$SIM_DIR"
mkdir -p "$BASE_PATH"

# Starting point determined by the second input
if [[ "$RESUME" == "y" ]]; then
    # Counting batches already in here
    LAST_BATCH=$(ls "$BASE_PATH" | grep "^batch_${SIM_DIR}_" | wc -l)
    START_BATCH=$LAST_BATCH
    echo "Starting from batch $START_BATCH"
else
    START_BATCH=0
    echo "Starting from scratch"
fi

# Launching the batches
for i in $(seq $START_BATCH $((START_BATCH + TOTAL_BATCHES - 1))); do
    BATCH_NAME="${SIM_DIR}/batch_${SIM_DIR}_${i}"
    OFFSET=$((i * SIMS_PER_BATCH))

    echo "Launching $BATCH_NAME with OFFSET=$OFFSET"

    sbatch \
      --export=ALL,BATCH_NAME=$BATCH_NAME,NUM_SIMS=$SIMS_PER_BATCH,OFFSET=$OFFSET \
      job_run_batch.slurm
done
