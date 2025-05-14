#!/bin/bash -l

# Echo des paramètres SLURM utiles (pour debug)
echo "BATCH=$BATCH_NAME, NUM_SIMS=$NUM_SIMS, OFFSET=$OFFSET"

#SBATCH --job-name=${BATCH_NAME}
#SBATCH -o $HOME/M2_Internship_CBGP/Project_Model/logs/%x_%j.out
#SBATCH -e $HOME/M2_Internship_CBGP/Project_Model/logs/%x_%j.err
#SBATCH -t 24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH -p long

# Activation de l'environnement conda
cd ~/M2_Internship_CBGP/Project_Model/Scripts/Python
source ~/miniconda3/etc/profile.d/conda.sh
conda activate stage_m2_env

# Lancement des simulations
for ((i=0; i<NUM_SIMS; i++)); do
    global_offset=$((OFFSET + i))
    echo "🚀 Simulation $global_offset in batch $BATCH_NAME"

    python main_cluster.py \
        --batch "$BATCH_NAME" \
        --num_simulations 1 \
        --name_prefix sim \
        --offset "$global_offset"
done
