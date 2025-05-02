#!/bin/bash -l
#SBATCH --job-name=sim_m2
#SBATCH --output=~/M2_Internship_CBGP/Project_Model/logs/sim_%j.out
#SBATCH --error=~/M2_Internship_CBGP/Project_Model/logs/sim_%j.err
#SBATCH --time=020000
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# Aller dans le dossier du dépôt GitHub
cd ~M2_Internship_CBGP/Project_Model/Scripts/Python

# Exécuter le script principal
python3 main_cluster.py --batch batch_test --num_simulations 5 --name_prefix sim
