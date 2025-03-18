# Packages needed to run python script
import os
import json
import time
import subprocess
import tskit                            # type: ignore
import pyslim                           # type: ignore
import msprime                          # type: ignore
import numpy as np                      # type: ignore
import matplotlib.pyplot as plt         # type: ignore
import warnings
import random

# Generating a new file for each simulation
timestamp = time.strftime("%Y%m%d_%H%M%S")
slurm_job_id = os.getenv("SLURM_JOB_ID", "local")
simulation_id = f"sim_{timestamp}_{slurm_job_id}"
output_folder = f"/scratch/$USER/slim_simulations/{simulation_id}/"
os.makedirs(output_folder, exist_ok = True)

# Creating JSON file with every simulation parameters
config = {
    "simulation_id" : simulation_id,
    "pop_size" : int(np.random.lognormal(mean=5, sigma=0.3)),
    "num_loci" : random.choice([5, 10, 20, 50]),
    "num_generations" : random.randint(5, 15),
    "sample_sizes" : [random.randint(20, 50), random.randint(20, 50)],
    "mutation_rate" : np.random.uniform(1e-5, 1e-2),
    "recap_Ne" : int(np.random.uniform(50, 500)),
    "output_folder" : output_folder,
    "timestamp" : timestamp,
    "seed" : random.randint(1, 10**6)
}

slim_config_file = os.path.join(output_folder, "slim_config.txt")
with open(slim_config_file, "w") as f:
    for key, value in config.items():
        if isinstance(value, list):
            value = ",".join(map(str, value))
        f.write(f"{key}={value}\n")

print(f"✅ Fichier de configuration SLiM généré : {slim_config_file}")

# Run SLiM script with generated parameters
slim_script = "Sim_model.slim"
slim_command = ["slim", "-d", f"config_file='{slim_config_file}'", slim_script]
with open(os.path.join(output_folder, "slim.log"), "w") as log_file:
    slim_process = subprocess.run(slim_command, stdout = log_file, stderr = log_file, text = True)
