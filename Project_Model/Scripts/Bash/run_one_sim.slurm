#!/bin/bash -l
#SBATCH --job-name=sim_test
#SBATCH --output=logs/sim_test_%j.out
#SBATCH --error=logs/sim_test_%j.err
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=long

# Aller dans le bon dossier
cd ~/M2_Internship_CBGP/Project_Model/Scripts/Python

# Activer l’environnement conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate stage_m2_env

# Lancer une seule simulation
python main_cluster.py \
    --batch test_debug \
    --num_simulations 1 \
    --name_prefix sim \
    --offset 0
