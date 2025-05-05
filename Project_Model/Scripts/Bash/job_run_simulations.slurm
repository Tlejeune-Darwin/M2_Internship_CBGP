#!/bin/bash -l
#SBATCH --job-name=${BATCH}
#SBATCH --output=sim_%A_%a.out
#SBATCH --error=sim_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --array=0-$((${SIMS}-1))%${MAX}

# Activer conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate stage_m2_env

# Aller dans le dossier du script
cd ~/M2_Internship_CBGP/Project_Model/Scripts

# Lancer la simulation
python main_cluster.py \
    --batch ${BATCH} \
    --num_simulations 1 \
    --name_prefix sim \
    --offset $((${OFFSET} + ${SLURM_ARRAY_TASK_ID}))

