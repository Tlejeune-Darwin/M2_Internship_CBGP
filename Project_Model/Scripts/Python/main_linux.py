from run_simulation_linux import run_simulation_linux  # Name of the simulation script
import argparse
import time


# --- Argparse module ---
# This part allows the simulation to be run n times
parser = argparse.ArgumentParser() # Convert to a command line
parser.add_argument("-n", "--num_simulations", type=int, default=10)
parser.add_argument("-d", "--directory", type=str, default="simulations", help="Top-level folder name to store all simulations.")
parser.add_argument("--pop_size", type=int, help="Census population size (optional).")
parser.add_argument("--num_loci", type=int, help="Number of loci (optional).")
parser.add_argument("--name_prefix", type=str, default="sim", help="Prefix for simulation folders (default: sim)")
args = parser.parse_args()

for i in range(args.num_simulations):  # 10 simulations
    run_simulation_linux(        
        base_dir=args.directory,
        pop_size=args.pop_size,
        num_loci=args.num_loci,
        sim_prefix=args.name_prefix
    )
    time.sleep(0.5)  # Latency between every simulation, necessary for RAM purpose
