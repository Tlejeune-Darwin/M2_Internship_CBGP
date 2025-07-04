from run_simulation_cluster import run_simulation_cluster
import argparse
import os

# --- Detects the user desktop directory (cross-platform) --- #
def get_desktop_path():
    desktop_fr = os.path.join(os.path.expanduser("~"), "Bureau") 
    desktop_en = os.path.join(os.path.expanduser("~"), "Desktop")
    return desktop_fr if os.path.isdir(desktop_fr) else (desktop_en if os.path.isdir(desktop_en) else os.path.expanduser("~"))

# --- Fallback prompt if arguments are missing --- #
# In cluster, this part is never showed because the bash script already got all the necessary informations. 
def ask_if_missing(value, label, type_func=str, default=None):
    if value is not None:
        return value
    msg = f"{label} [{'default: ' + str(default) if default is not None else 'required'}]: "
    inp = input(msg)
    return type_func(inp) if inp else default

# --- Command-line arguments --- #
# Command lines obtained from the bash script
# If not, refer to the names in each argument
parser = argparse.ArgumentParser(description="Launch a batch of SLiM simulations.")
parser.add_argument("--batch", type=str, help="Name of the subfolder grouping the batch of simulations.")
parser.add_argument("-n", "--num_simulations", type=int, help="Number of simulations to run.")
parser.add_argument("--offset", type=int, default=0, help="Starting index offset for sim numbering")
parser.add_argument("--name_prefix", type=str, default="sim", help="Prefix for simulation folder names (default: sim).")
args = parser.parse_args()

# --- Fill missing arguments interactively --- #
args.batch = ask_if_missing(args.batch, "Batch name (simulation folder name)", str)
args.num_simulations = ask_if_missing(args.num_simulations, "Number of simulations", int, 1)

# --- Create the base simulation folder on the desktop --- #
# In cluster, this part is put inside the "results" folder
global_dir, batch_dir = args.batch.split("/", 1)
base_results_dir = os.path.expanduser(f"~/results/{global_dir}/{batch_dir}")
os.makedirs(base_results_dir, exist_ok=True)
sim_base_dir = base_results_dir

# --- Run the simulations --- #
for i in range(args.num_simulations):
    run_simulation_cluster(
        base_dir=sim_base_dir,  # Pass the absolute path
        sim_prefix=args.name_prefix,
        offset=args.offset + i
    )

