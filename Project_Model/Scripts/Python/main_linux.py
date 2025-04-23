from run_simulation_linux import run_simulation_linux  # Name of the simulation script
import argparse
import time
import os

def ask_if_missing(value, label, type_func=str, default=None):
    if value is not None:
        return value
    msg = f"{label} [{'default: ' + str(default) if default else 'required'}]: "
    inp = input(msg)
    return type_func(inp) if inp else default

parser = argparse.ArgumentParser()
parser.add_argument("--batch", type=str, help="Nom du sous-dossier regroupant un ensemble de simulations.")
parser.add_argument("-n", "--num_simulations", type=int, help="Nombre de simulations à lancer")
parser.add_argument("--pop_size", type=int)
parser.add_argument("--num_loci", type=int)

args = parser.parse_args()

# Demande interactive si les arguments sont absents
args.batch = ask_if_missing(args.batch, "Nom du batch (dossier de simulation)")
args.num_simulations = ask_if_missing(args.num_simulations, "Nombre de simulations", int, 1)
args.pop_size = ask_if_missing(args.pop_size, "Taille de population (laisser vide pour aléatoire)", int, None)
args.num_loci = ask_if_missing(args.num_loci, "Nombre de loci", int, 20)

# Exécution
for i in range(args.num_simulations):
    run_simulation_linux(
        base_dir=os.path.join("simulations", args.batch),
        pop_size=args.pop_size,
        num_loci=args.num_loci,
        sim_prefix="sim"
    )

# --- Argparse module ---
# This part allows the simulation to be run n times
parser = argparse.ArgumentParser() # Convert to a command line
parser.add_argument("-n", "--num_simulations", type=int, default=10)
parser.add_argument("--batch", type=str, required=True,
                    help="Nom du sous-dossier regroupant un ensemble de simulations.")
args = parser.parse_args()

for i in range(args.num_simulations):  # 10 simulations
    run_simulation_linux(        
        base_dir=os.path.join("simulations", args.batch),
        pop_size=args.pop_size,
        num_loci=args.num_loci,
        sim_prefix=args.name_prefix
    )
    time.sleep(0.001)  # Latency between every simulation, necessary for RAM purpose
