#!/bin/bash -l
#SBATCH --job-name=sim_array
#SBATCH --output=sim_%A_%a.out
#SBATCH --error=sim_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=long
#SBATCH --array=0-99%4

# Activer conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate stage_m2_env

# Aller dans le dossier du script
cd ~/M2_Internship_CBGP/Project_Model/Scripts/Python

# Lancer 1 simulation avec l'indice du tableau comme offset
python main_cluster.py --batch batch_array --num_simulations 1 --name_prefix sim --offset ${SLURM_ARRAY_TASK_ID}
