import tskit # type: ignore
import pyslim # type: ignore
import msprime # type: ignore
import warnings # type: ignore
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
import IPython.display # type: ignore
import os
import subprocess
import shutil
import sys
import os

                                                    # ---___---___---___--- Directory ---___---___---___--- #

# Define a base file where everything will be run
BASE_DIR = os.path.dirname(os.path.abspath(__file__))  # Path of the currently running Python script
SLIM_SCRIPT = os.path.join(BASE_DIR, "Model_nWF_10.slim")
OUTPUT_DIR = os.path.join(BASE_DIR, "output_trees")
TREE_FILE = os.path.join(OUTPUT_DIR, "simulation.trees")
RECAP_FILE = os.path.join(OUTPUT_DIR, "simulation_recap.trees")

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"Base directory : {BASE_DIR}") # PROMPT : can be ignored
print(f"SLiM script : {SLIM_SCRIPT}") # PROMPT : can be ignored
print(f"Output directory : {OUTPUT_DIR}") # PROMPT : can be ignored

# Define log output path
log_file_path = os.path.join(OUTPUT_DIR, "simulation_output.log")

# Write the output inside de log file
log_file = open(log_file_path, "w")

# Commands and prompts in console will be written in the log file
sys.stdout = log_file

                                                    # ---___---___---___--- Parameters ---___---___---___--- #

RECAP_Ne = 200  # Effective population size for recapitation
RECOMB_RATE = 0.5  # Recombination rate
MUT_RATE = 1e-2  # Mutation rate
HI = 50
LO = 5

                                                    # ---___---___---___--- Run SLiM ---___---___---___--- #

# Automatically find SLiM on the computer
slim_path = shutil.which("slim")
if slim_path is None:
    raise FileNotFoundError("SLiM not found. Please define SLiM_PATH.") # PROMPT : can be ignored

# Make the dynamic command
log_file = os.path.join(OUTPUT_DIR, "simulation.log")
command = f'"{slim_path}" "{SLIM_SCRIPT}" > "{log_file}" 2>&1'

print(f"Launching SLiM : {command}") # PROMPT : can be ignored

if not os.path.exists(SLIM_SCRIPT):
    raise FileNotFoundError(f"Error : SLiM script {SLIM_SCRIPT} doesn't exist !") # PROMPT : can be ignored
else:
    print(f"SLiM script detected : {SLIM_SCRIPT}") # PROMPT : can be ignored

def run_slim():
    """ Run SLiM with generic path """
    log_file = os.path.join(OUTPUT_DIR, "simulation.log")
    command = f'"{slim_path}" "{SLIM_SCRIPT}" > "{log_file}" 2>&1'

    print(f"Run SLiM simulation : {command}") # PROMPT : can be ignored 

    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Simulation over, log saved to {log_file}") # PROMPT : can be ignored
    except subprocess.CalledProcessError as e:
        print(f"Error while executing SLiM : {e}") # PROMPT : can be ignored

if __name__ == "__main__":
    run_slim()

# Delete the TimeUnitMismatch Warning 
warnings.simplefilter("ignore", msprime.TimeUnitsMismatchWarning)

                                                    # ---___---___---___--- Filtering ---___---___---___--- #

# Load the tree sequence made by SLiM
tree_sequence = tskit.load("output_trees/simulation.trees")

# Identify "Remember" and "Retained" individuals
kept_individuals = []
for ind in tree_sequence.individuals():
    if (ind.flags & pyslim.INDIVIDUAL_REMEMBERED) or (ind.flags & pyslim.INDIVIDUAL_RETAINED):
        kept_individuals.append(ind.id)

print(f"Number of Remembered and Retained individuals {len(kept_individuals)}") # PROMPT : can be ignored

# Take the nodes from the filtered individuals only
kept_nodes = []
for ind in kept_individuals:
    kept_nodes.extend(tree_sequence.individual(ind).nodes)

print(f"Total number of nodes after filtering {len(kept_nodes)}") # PROMPT : can be ignored

# Simplify tree sequence to keep the filtered individuals
filtered_ts = tree_sequence.simplify(kept_nodes, keep_input_roots=True)

print(f"Number of trees after filtering {filtered_ts.num_trees}") # PROMPT : can be ignored

# Verify simplification worked
if filtered_ts.num_trees == 0:
    print("No tree found after simplification. Please verify individual flags.") # PROMPT : can be ignored

                                                    # ---___---___---___--- Recapitation ---___---___---___--- #

# Define demography
demography = msprime.Demography()
demography.add_population(name="p1", initial_size=RECAP_Ne)

# Recapitate the simplified tree
recap_ts = pyslim.recapitate(filtered_ts, recombination_rate=RECOMB_RATE, demography=demography)

# New simplification to avoid more than one root
recap_ts = recap_ts.simplify(keep_input_roots=False)

# Verify that recapitation worked
for tree in recap_ts.trees():
    print(f"After recapitation - Root number : {tree.num_roots}") # PROMPT : can be ignored

# Text tree in console for verification
##         print("\nTree plot - Text version :\n") # PROMPT : can be ignored          ##
##         print(recap_ts.first().draw_text())         ##

# Most Recent Common Ancestor (MRCA)
mrca_ages = []
for tree in recap_ts.trees():
    mrca_age = tree.time(tree.root)  # Most recent common ancestor age
    mrca_ages.append(mrca_age)

for i, tree in enumerate(recap_ts.trees()):
    print(f"Tree {i+1} : MRCA Age = {tree.time(tree.root):.2f} generations")


                                                    # ---___---___---___--- Add Mutation ---___---___---___--- #

# Adding mutation after recapitation
rates = np.linspace(1e-2, 1e-4, num=100)
variances = []
mut_model = msprime.SMM(lo=LO, hi=HI)  # Stepwise Mutation Model with shifted values
mut_ts = msprime.sim_mutations(recap_ts, rate=MUT_RATE, model=mut_model, random_seed=42)

# Function that convert the output into some more useful data
def copy_number_matrix(ts):
    """
    Returns the copy number matrix from the specified tree sequence
    simulated under a MicrosatMutationModel.
    """
    C = np.zeros((ts.num_sites, ts.num_samples), dtype=int)  # Copy number matrix
    
    for var in ts.variants():
        print(f"Site {var.site.id}: {var.alleles}")  # Verify allele list before conversion PROMPT : can be ignored
        # Convert allele in integer
        alleles = np.array([
            int(allele) if allele is not None else 0  # Replace "None" by "0", if not stop the script
            for allele in var.alleles
        ])
        
        C[var.site.id] = alleles[var.genotypes]  # Assign alleles to genotypes
    
    return C

# Graph showing variance in the number of repeat with the variation in mutation rate
for r in rates:
    mts = msprime.sim_mutations(recap_ts, rate=r, model=mut_model, random_seed=1)
    print(f"Mutation rate : {r}, Total number of mutations : {mts.num_mutations}") # PROMPT : can be ignored
    C = copy_number_matrix(mts)
    print(f"Sample Matrix C :\n", C[:5])  # Show the first five lines - PROMPT : can be ignored
    if C.size == 0 or np.all(C == 0):  # Verify if matrix is empty or full of 0
        print(f"No mutations, r={r}. Ignored.") #PROMPT : can be ignored
        variances.append(np.nan)  # Add NaN values if necessary to avoid errors
    else:
        variances.append(C.var())

window_size = 5
smoothed_variances = np.convolve(variances, np.ones(window_size)/window_size, mode='valid')

plt.plot(rates, variances, alpha=0.3, label="Raw data")
plt.plot(rates[window_size-1:], smoothed_variances, color="red", label="Smoothed (Moving Avg)")
plt.xlabel('Mutation rate')
plt.ylabel('Variance in repeat number')
plt.title("Mutation rate impact on repeat variance")
plt.show()

print(f"Total number of added mutations : {mut_ts.num_mutations}") # PROMPT : can be ignored

mutation_counts = [len(site.mutations) for site in mut_ts.sites()]
unique_mut_counts, site_frequencies = np.unique(mutation_counts, return_counts=True)

plt.bar(unique_mut_counts, site_frequencies, edgecolor="black", alpha=0.7)
plt.xlabel("Number of mutations per site")  # Axe X = Mutation number
plt.ylabel("Number of sites")  # Axe Y = Site that have this number of mutation
plt.title("Mutation number distribution per site")
plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.show()

# Verify first tree content
tree = mut_ts.first()
print(f"Node number in the final tree : {len(tree.preorder())}") # PROMPT : can be ignored
print(f"Sample number after simplification {tree.num_samples()}") # PROMPT : can be ignored

# Extraction of copy numbers
copy_numbers = copy_number_matrix(mut_ts)

# First line display after verification
print("Extract of copy number matrix :") # PROMPT : can be ignored
print(copy_numbers[:5])  # Display 5 first sites - PROMPT : can be ignored
print(f"Shape of copy_numbers matrix: {copy_numbers.shape}")  # PROMPT : can be ignored

output_gen_file = os.path.join(OUTPUT_DIR, "simulation_data.gen")

# Get individual number and loci
copy_numbers = copy_numbers.T # T to transpose line and column 
num_individuals, num_loci = copy_numbers.shape 
loci_names = " ".join([f"Locus_{i+1}" for i in range(num_loci)]) # Create the list of loci in a specific combination

# To make Neestimator work, it is required to write the log file in a very specific way (For GENEPOP format, it is different in Fstat format) 
with open(output_gen_file, "w") as f:
    f.write("Simulated GENEPOP Data \n")  # Title
    f.write(loci_names + "\n") # Loci number
    f.write("pop\n") # Population name

    for i in range(num_individuals):
        genotype_line = ", " + " ".join(f"{int(copy_numbers[i, j]):03}{int(copy_numbers[i, j]):03}" for j in range(num_loci))  # Convert the alleles and double it (diploid) to get the genotype of each individual
        f.write(f"Indiv_{i+1} {genotype_line}\n")  # Individual ID + their genotype

print(f".gen file generated : {output_gen_file}") #PROMPT : can be ignored

                                                    # ---___---___---___--- Tree Display ---___---___---___--- #

# Generat SVG
svg_output = mut_ts.draw_svg(
    size=(3000, 3000),
    node_labels={},
    mutation_labels=None,
    time_scale="log_time",  # Using "Rank" is better to visualize all ancestor by recapitation and mutations but less representative
    x_axis=True, # Genomic position
    y_axis=True  # Time axis
)

# Graphical representation of mutation rate
from matplotlib import pyplot as plt # type: ignore

# Saved in output_trees file
tree_file_path = os.path.join(OUTPUT_DIR, "tree_output.svg")
with open(tree_file_path, "w") as f:
    f.write(svg_output)

print(f"Tree recapitated and saved to {tree_file_path}")  # PROMPT : can be ignored

