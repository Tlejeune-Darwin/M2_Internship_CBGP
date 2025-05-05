#!/bin/bash -l

# --- Sécurités ---
: "${OFFSET:?OFFSET non défini}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID non défini}"
: "${BATCH:?Nom du batch manquant}"

#SBATCH --job-name=${BATCH}
#SBATCH --output=sim_%A_%a.out
#SBATCH --error=sim_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# --- Active conda (modifie ce chemin si besoin) ---
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate stage_m2_env

# --- Aller dans le bon répertoire ---
cd "$HOME/M2_Internship_CBGP/Project_Model/Scripts/Python"

# --- Lancer le script ---
python main_cluster.py \
  --batch "${BATCH}" \
  --num_simulations 1 \
  --name_prefix sim \
  --offset $((OFFSET + SLURM_ARRAY_TASK_ID))
