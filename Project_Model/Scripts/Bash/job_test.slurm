#!/bin/bash -l
#SBATCH --job-name=test_env
#SBATCH --output=test_env_%j.out
#SBATCH --error=test_env_%j.err
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --partition=short

# Activer conda proprement
source ~/miniconda3/etc/profile.d/conda.sh
conda activate stage_m2_env

# Afficher infos utiles
echo "Dossier : $(pwd)"
echo "Python utilisé : $(which python)"
python --version
