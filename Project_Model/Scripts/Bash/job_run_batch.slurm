#!/bin/bash -l

# Calcul manuel de l’intervalle d’array
ARRAY_MAX=$(($SIMS - 1))

#SBATCH --job-name=batch_run
#SBATCH --output=sim_%A_%a.out
#SBATCH --error=sim_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --array=0-9999%10

cd $HOME/M2_Internship_CBGP/Project_Model/Scripts/Bash/

# Lancer une simulation unique avec offset correct
python main_cluster.py \
  --batch ${BATCH} \
  --num_simulations 1 \
  --name_prefix sim \
  --offset $((${OFFSET:=0} + ${SLURM_ARRAY_TASK_ID}))
