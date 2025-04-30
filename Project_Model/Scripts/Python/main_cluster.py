from run_simulation_cluster import run_simulation_cluster
import argparse
import os

# --- Detects the user desktop directory (cross-platform) --- #
def get_desktop_path():
    desktop_fr = os.path.join(os.path.expanduser("~"), "Bureau")
    desktop_en = os.path.join(os.path.expanduser("~"), "Desktop")
    return desktop_fr if os.path.isdir(desktop_fr) else (desktop_en if os.path.isdir(desktop_en) else os.path.expanduser("~"))

# --- Fallback prompt if arguments are missing --- #
def ask_if_missing(value, label, type_func=str, default=None):
    if value is not None:
        return value
    msg = f"{label} [{'default: ' + str(default) if default is not None else 'required'}]: "
    inp = input(msg)
    return type_func(inp) if inp else default

# --- Command-line arguments --- #
parser = argparse.ArgumentParser(description="Launch a batch of SLiM simulations.")
parser.add_argument("--batch", type=str, help="Name of the subfolder grouping the batch of simulations.")
parser.add_argument("-n", "--num_simulations", type=int, help="Number of simulations to run.")
parser.add_argument("--name_prefix", type=str, default="sim", help="Prefix for simulation folder names (default: sim).")
args = parser.parse_args()

# --- Fill missing arguments interactively --- #
args.batch = ask_if_missing(args.batch, "Batch name (simulation folder name)", str)
args.num_simulations = ask_if_missing(args.num_simulations, "Number of simulations", int, 1)

# --- Create the base simulation folder on the desktop --- #
base_results_dir = os.path.expanduser("results")
sim_base_dir = os.path.join(base_results_dir, args.batch)
os.makedirs(sim_base_dir, exist_ok=True)

# --- Run the simulations --- #
for i in range(args.num_simulations):
    run_simulation_cluster(
        base_dir=sim_base_dir,  # Pass the absolute path
        sim_prefix=args.name_prefix
    )

