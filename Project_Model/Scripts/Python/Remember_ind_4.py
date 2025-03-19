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

RECAP_Ne = 6  # Effective population size for recapitation
RECOMB_RATE = 0.5  # Recombination rate
MUT_RATE = 1e-1  # Mutation rate
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

for ind in filtered_ts.individuals():
    print(f"Individu {ind.id} -> Nœuds {ind.nodes}")



print(f"Number of trees after filtering {filtered_ts.num_trees}") # PROMPT : can be ignored

# Verify simplification worked
if filtered_ts.num_trees == 0:
    print("No tree found after simplification. Please verify individual flags.") # PROMPT : can be ignored

# Nombre total d'individus dans l'arbre initial
total_individuals = tree_sequence.num_individuals
print(f"Nombre total d'individus avant filtrage : {total_individuals}")

# Nombre d'individus retenus
print(f"Nombre d'individus conservés : {len(kept_individuals)}")

# Nombre total de nœuds dans l'arbre initial
total_nodes = tree_sequence.num_nodes
print(f"Nombre total de nœuds avant filtrage : {total_nodes}")

# Nombre de nœuds conservés
print(f"Nombre de nœuds conservés après filtrage : {len(kept_nodes)}")


print(f"ID des individus retenus : {kept_individuals}")
for ind_id in kept_individuals:
    print(f"Individu {ind_id}, nœuds associés : {tree_sequence.individual(ind_id).nodes}")


valid_nodes = {node for ind in kept_individuals for node in tree_sequence.individual(ind).nodes}
node_check = all(node in valid_nodes for node in kept_nodes)

print(f"Tous les nœuds retenus appartiennent bien à des individus filtrés : {node_check}")



# Affichage d'un arbre avant filtrage
output_tree_file = os.path.join(OUTPUT_DIR, "tree_before_filtering.txt")
with open(output_tree_file, "w", encoding="utf-8") as f:
    f.write(tree_sequence.first().draw_text())

print(f"L'arbre avant filtrage a été enregistré dans {output_tree_file}")


# Affichage d'un arbre après filtrage
output_filtered_tree_file = os.path.join(OUTPUT_DIR, "tree_after_filtering.txt")
with open(output_filtered_tree_file, "w", encoding="utf-8") as f:
    f.write(filtered_ts.first().draw_text())

print(f"L'arbre après filtrage a été enregistré dans {output_filtered_tree_file}")



remaining_inds = set(ind.id for ind in filtered_ts.individuals())
missing_inds = set(kept_individuals) - remaining_inds
print(f"Individus qui ont été supprimés après filtrage : {missing_inds}")

for tree in filtered_ts.trees():
    if tree.num_roots > 1:
        print(f"Avertissement : Arbre avec plusieurs racines détecté après simplification ({tree.num_roots} racines).")




                                                    # ---___---___---___--- Recapitation ---___---___---___--- #

# Define demography
demography = msprime.Demography()
demography.add_population(name="p1", initial_size=RECAP_Ne)

# Recapitate the simplified tree
recap_ts = pyslim.recapitate(filtered_ts, recombination_rate=RECOMB_RATE, demography=demography)

# New simplification to avoid more than one root
#recap_ts = recap_ts.simplify(keep_input_roots=False)

# Verify that recapitation worked
for tree in recap_ts.trees():
    print(f"After recapitation - Root number : {tree.num_roots}") # PROMPT : can be ignored

# Text tree in console for verification
print("\nTree plot - Text version :\n") # PROMPT : can be ignored          ##
output_tree_file = os.path.join(OUTPUT_DIR, "recapitated_tree.txt")
with open(output_tree_file, "w", encoding="utf-8") as f:
    f.write(recap_ts.first().draw_text())

print(f"L'arbre recapité a été enregistré dans {output_tree_file}. Ouvre-le avec un éditeur UTF-8.")
        

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

    for i in range(0, num_individuals, 2):  # On saute une ligne sur deux
        genotype_line = ", " + " ".join(f"{int(copy_numbers[i, j]):03}{int(copy_numbers[i+1, j]):03}" for j in range(num_loci))
        f.write(f"Indiv_{(i//2)+1}{genotype_line}\n")  # Unique ID per individual

for mut in mut_ts.mutations():
    node = mut_ts.node(mut.node)
    print(f"Mutation {mut.derived_state} attachée au nœud {mut.node} (temps {node.time})")


print(f".gen file generated : {output_gen_file}") #PROMPT : can be ignored


                                                    # ---___---___---___--- Tree Display ---___---___---___--- #

def get_node_labels(ts, site_index=0):
    """
    Retourne un dictionnaire avec les labels des nœuds,
    contenant l'ID de l'individu et son génotype pour un site donné.
    """
    node_labels = {}

    # Obtenir la matrice de génotypes
    copy_numbers = copy_number_matrix(ts)

    # Associer chaque individu à son/ses nœuds et à son génotype
    for ind in ts.individuals():
        ind_id = ind.id  # ID de l'individu
        nodes = ind.nodes  # Noeuds associés

        # Vérifier qu’on a bien les informations génétiques
        if ind_id < copy_numbers.shape[1]:  
            for i, node in enumerate(nodes):  # Chaque nœud a son propre allèle
                genotype = copy_numbers[site_index, node]  # Génotype spécifique au nœud
                node_labels[node] = f"N:{node} | ID:{ind_id} | G:{genotype}"  # Afficher seulement cet allèle

    # Ajouter les labels pour les nœuds intermédiaires (ancêtres communs)
    for tree in ts.trees():
        for node in tree.nodes():
            if node not in node_labels:  # Ajouter seulement si ce n'est pas une feuille
                node_labels[node] = f"N:{node}"  # Afficher uniquement l’ID du nœud

    return node_labels


# Générer les labels pour l'affichage
node_labels = get_node_labels(mut_ts, site_index=0)  # Étudie le site 0 par défaut

def get_mutation_labels(ts):
    """
    Retourne un dictionnaire contenant les labels des mutations, indiquant si elles 
    représentent un gain (+1) ou une perte (-1), même si le nœud parent est un nœud recapité.
    """
    mutation_labels = {}

    for tree in ts.trees():
        last_known_state = {}  # Dictionnaire stockant la dernière mutation connue pour chaque lignée

        for mut in tree.mutations():
            node = mut.node  # Nœud portant la mutation
            current_allele = int(mut.derived_state)  # État après mutation

            # Trouver la mutation connue la plus proche dans l'ascendance
            previous_allele = None
            parent_node = node  # On commence par le parent immédiat
            
            while parent_node != -1:
                # Vérifier si une mutation précédente existe sur cette lignée
                if parent_node in last_known_state:
                    previous_allele = last_known_state[parent_node]
                    break  # On a trouvé un état génétique précédent
                
                parent_node = tree.parent(parent_node)  # Remonter encore si nécessaire
            
            # Si aucune mutation connue n'a été trouvée, ignorer la mutation (évite les erreurs)
            if previous_allele is not None:
                mutation_effect = current_allele - previous_allele
                sign = "+" if mutation_effect > 0 else "-"

                mutation_labels[mut.id] = f"{sign}{abs(mutation_effect)}"  # Ex : "+1" ou "-1"

            # Stocker la mutation actuelle pour qu'elle serve de référence aux suivantes
            last_known_state[node] = current_allele

    return mutation_labels




mutation_labels = get_mutation_labels(mut_ts)

# Génération du SVG avec les labels
svg_output = mut_ts.draw_svg(
    size=(3000, 1000),
    node_labels=node_labels,  # On ajoute les labels
    mutation_labels=mutation_labels,
    time_scale="rank",
    x_axis=True,
    y_axis=True
)

# Sauvegarde de l'arbre annoté
tree_file_path = os.path.join(OUTPUT_DIR, "tree_with_labels.svg")
with open(tree_file_path, "w") as f:
    f.write(svg_output)

print(f"Tree with labels saved to {tree_file_path}")


