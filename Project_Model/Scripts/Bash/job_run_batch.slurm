#!/bin/bash -l

#SBATCH -J Job_simu
#SBATCH -p long
#SBATCH -N 2
#SBATCH -n 32
#SBATCH -A stage_m2
#SBATCH -t 2-12
#SBATCH -o ./%N.out
#SBATCH -e ./%N.err

# Conda environment activation
cd ~/M2_Internship_CBGP/Project_Model/Scripts/Python
source ~/miniconda3/etc/profile.d/conda.sh
conda activate stage_m2_env

# Launching simulation (from the python script)
for ((i=0; i<NUM_SIMS; i++)); do
    global_offset=$((OFFSET + i))
    echo "🚀 Simulation $global_offset in batch $BATCH_NAME"

    python main_cluster.py \
        --batch "$BATCH_NAME" \
        --num_simulations 1 \
        --name_prefix sim \
        --offset "$global_offset"
done
