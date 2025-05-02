#!/bin/bash -l

# Sécurités : arrêt immédiat si une variable est absente
: "${OFFSET:?OFFSET non défini}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID non défini}"
: "${BATCH:?Nom du batch manquant}"

# Job SLURM
#SBATCH --job-name=${BATCH}
#SBATCH --output=sim_%A_%a.out
#SBATCH --error=sim_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --array=0-9%10

# Aller dans le bon dossier
cd "$HOME/M2_Internship_CBGP/Project_Model/Scripts/Python"

# Lancer la simulation
python main_cluster.py \
  --batch "${BATCH}" \
  --num_simulations 1 \
  --name_prefix sim \
  --offset $((${OFFSET} + ${SLURM_ARRAY_TASK_ID}))
