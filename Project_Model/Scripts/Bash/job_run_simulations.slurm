#!binbash
#SBATCH --job-name=sim_m2
#SBATCH --output=logssim_%j.out
#SBATCH --error=logssim_%j.err
#SBATCH --time=020000
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# Optionnel  purge et chargement de modules si nécessaire
module purge

# Aller dans le dossier du dépôt GitHub
cd ~ScriptPython

# Exécuter le script principal
python main_cluster.py --batch batch_test --num_simulations 5 --name_prefix sim
