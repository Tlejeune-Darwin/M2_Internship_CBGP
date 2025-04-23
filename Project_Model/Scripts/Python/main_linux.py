from run_simulation_linux import run_simulation_linux  # Name of the simulation script
import argparse
import time
import os

def ask_if_missing(value, label, type_func=str, default=None):
    if value is not None:
        return value
    msg = f"{label} [{'default: ' + str(default) if default is not None else 'required'}]: "
    inp = input(msg)
    return type_func(inp) if inp else default

# --- Arguments de ligne de commande --- #
parser = argparse.ArgumentParser(description="Lance un batch de simulations SLiM.")
parser.add_argument("--batch", type=str, help="Nom du sous-dossier regroupant un ensemble de simulations.")
parser.add_argument("-n", "--num_simulations", type=int, help="Nombre de simulations à lancer")
parser.add_argument("--name_prefix", type=str, default="sim", help="Préfixe pour le nom des simulations (default: sim)")
args = parser.parse_args()

# --- Complétion interactive si nécessaire --- #
args.batch = ask_if_missing(args.batch, "Nom du batch (dossier de simulation)", str)
args.num_simulations = ask_if_missing(args.num_simulations, "Nombre de simulations", int, 1)

# --- Exécution du batch de simulations --- #
for i in range(args.num_simulations):
    run_simulation_linux(
        base_dir=os.path.join("simulations", args.batch),
        sim_prefix=args.name_prefix
    )
    time.sleep(0.001)  # Pause minimale pour éviter collisions et surcharge CPU
