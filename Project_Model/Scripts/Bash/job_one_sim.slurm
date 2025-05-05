#!/bin/bash -l
#SBATCH --job-name=sim_m2
#SBATCH --output=sim_m2_%j.out
#SBATCH --error=sim_m2_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --partition=long

cd ~/M2_Internship_CBGP/Project_Model/Scripts/Python
source ~/miniconda3/etc/profile.d/conda.sh
conda activate stage_m2_env

python main_cluster.py --batch test_debug --num_simulations 1 --name_prefix sim
