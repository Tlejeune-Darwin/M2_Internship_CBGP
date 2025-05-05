#!/bin/bash -l

echo "BATCH=$BATCH_NAME, NUM_SIMS=$NUM_SIMS, OFFSET=$OFFSET"

#SBATCH --job-name=${BATCH_NAME}
#SBATCH --output=$HOME/GitHub/M2_Internship_CBGP/Project_Model/logs/%x_%j.out
#SBATCH --error=$HOME/GitHub/M2_Internship_CBGP/Project_Model/logs/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=long

cd ~/M2_Internship_CBGP/Project_Model/Scripts/Python
source ~/miniconda3/etc/profile.d/conda.sh
conda activate stage_m2_env

# Simulation loop
for i in $(seq 0 $((NUM_SIMS - 1))); do
    global_offset=$((OFFSET + i))
    echo "Simulation $global_offset dans $BATCH_NAME"
    python main_cluster.py \
        --batch "$BATCH_NAME" \
        --num_simulations 1 \
        --name_prefix sim \
        --offset "$global_offset"
done
