#!/bin/bash -l

#SBATCH --job-name=batch_100
#SBATCH --output=sim_100_%j.out
#SBATCH --error=sim_100_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=long

# Aller dans le dossier du script
cd ~/M2_Internship_CBGP/Project_Model/Scripts/Python

# Activer conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate stage_m2_env

# Paramètres du batch
BATCH_NAME="batch_100"
NUM_SIMS=100
NAME_PREFIX="sim"

# Boucle sur les simulations
for i in $(seq 0 $((NUM_SIMS - 1))); do
    echo "Lancement de la simulation $i"
    python main_cluster.py \
        --batch "$BATCH_NAME" \
        --num_simulations 1 \
        --name_prefix "$NAME_PREFIX" \
        --offset "$i"
done
