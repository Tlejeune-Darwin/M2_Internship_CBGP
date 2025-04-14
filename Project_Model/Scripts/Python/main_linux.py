from run_simulation_linux import run_simulation_linux  # Name of the simulation script
import argparse
import time

# --- Argparse module ---
# This part allows the simulation to be run n times
parser = argparse.ArgumentParser() # Convert to a command line
parser.add_argument(  # Create the argument that will be converted
    "-n", "--num_simulations",
    type=int, 
    default=10,  # Number of simulations if nothing is specified
)
args = parser.parse_args()

for i in range(args.num_simulations):  # 10 simulations
    run_simulation_linux()  # Python script 
    time.sleep(0.5)  # Latency between every simulation, necessary for RAM purpose
