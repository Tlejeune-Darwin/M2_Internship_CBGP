#!/bin/bash -l
#SBATCH --job-name=sim_m2
#SBATCH --output=~sim_%j.out
#SBATCH --error=~sim_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=long

# Aller dans le dossier contenant le script
cd ~/M2_Internship_CBGP/Project_Model/Scripts

# Activer conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate stage_m2_env

# Exécuter ton script
python main_cluster.py --batch test --num_simulations 3 --name_prefix sim
