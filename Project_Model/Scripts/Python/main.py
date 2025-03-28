import os
from run_simulation import run_simulation  # Tu dois avoir une fonction principale dans run_simulation.py
import time
import argparse

# --- Lecture du nombre de simulations en argument ---
parser = argparse.ArgumentParser(description="Lancer des simulations évolutives")
parser.add_argument(
    "-n", "--num_simulations",
    type=int,
    default=10,
    help="Nombre de simulations à exécuter (défaut: 10)"
)
args = parser.parse_args()

def generate_simulation_parameters(n):
    """
    Génère des paramètres de simulation aléatoires ou contrôlés.
    Tu peux remplacer ça par la lecture d’un fichier JSON ou autre.
    """
    from random import randint, uniform
    return {
        "pop_size": randint(100, 500),
        "num_loci": 10,
        "num_generations": 20,
        "sample_size_1": "20",
        "sample_size_2": "30",
        "mutation_rate": round(uniform(0.001, 0.01), 4)
    }

for i in range(args.num_simulations):  # 10 simulations
    sim_params = generate_simulation_parameters(i)
    print(f"\n🚀 Simulation {i+1} — paramètres : {sim_params}")
    run_simulation(sim_params)  # Tu dois adapter cette fonction
    time.sleep(1)  # Laisse respirer un peu le système
