#!/bin/bash

# Nombre de batchs de 10 000 simulations chacun
NB_BATCHES=2
SIMS_PER_BATCH=10
MAX_CONCURRENT=10   # Nombre max de jobs en simultané

# Dossier contenant le script SLURM
SCRIPT_DIR="$HOME/M2_Internship_CBGP/Project_Model/Scripts/Bash"

for i in $(seq -w 0 $((NB_BATCHES - 1))); do
    BATCH_NAME="batch_$i"
    OFFSET=$((10#$i * SIMS_PER_BATCH))
    
    sbatch --export=BATCH=$BATCH_NAME,OFFSET=$OFFSET,SIMS=$SIMS_PER_BATCH,MAX=$MAX_CONCURRENT \
           "$SCRIPT_DIR/job_run_batch.slurm"
done