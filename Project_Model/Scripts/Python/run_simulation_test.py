def run_simulation_test():

    # Packages needed to run python script
    import os
    import subprocess
    import tskit                            # type: ignore
    import pyslim                           # type: ignore
    import msprime                          # type: ignore
    import numpy as np                      # type: ignore
    import matplotlib.pyplot as plt         # type: ignore
    import random
    import warnings
    import shutil
    import pandas as pd                     # type: ignore
    import re
    from datetime import datetime

    # ---___---___---___--- Setup directories ---___---___---___--- #

    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

    # Simulation folder
    all_simulations = os.path.join(SCRIPT_DIR, "simulations")
    os.makedirs(all_simulations, exist_ok=True)

    # Create a folder for each simulation
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    sim_id = f"sim_{timestamp}_local" # Give a number to the simulation's name
    sim_folder = os.path.join(all_simulations, sim_id)
    os.makedirs(sim_folder, exist_ok=True)

    # Creating text file with every simulation parameters
    config = {
        "simulation_id" : sim_id,
        "pop_size" : 200,
        "num_loci" : 19,
        "num_generations" : 9,
        "sample_sizes_Ne" : [50, 50],
        "sample_sizes_CMR" : [50, 50],
        "mutation_rate" : 5*10e-4,
        "recap_Ne" : 200,
        "output_folder" : sim_folder,
        "timestamp" : timestamp,
        "seed" : random.randint(1, 10**6)
    }

    slim_config_file = os.path.join(sim_folder, "slim_config.txt")
    with open(slim_config_file, "w") as f:
        for key, value in config.items():
            if isinstance(value, list):
                value = ",".join(map(str, value))
            f.write(f"{key}={value}\n")

    # ---___---___---___--- Run SLiM ---___---___---___--- #

    slim_script = os.path.join(SCRIPT_DIR, "Sim_model.slim")
    slim_executable = os.path.join(SCRIPT_DIR, "slim.exe")
    slim_config_file = slim_config_file.replace("\\", "/")
    slim_command = [slim_executable, "-d", f'config_file="{slim_config_file}"', slim_script]

    log_file_path = os.path.join(sim_folder, "slim.log")

    with open(log_file_path, "w") as log_file:
        slim_process = subprocess.run(slim_command, stdout = log_file, stderr = log_file, text = True, cwd = SCRIPT_DIR)

    # Delete the TimeUnitMismatch Warning 
    warnings.simplefilter("ignore", msprime.TimeUnitsMismatchWarning)

    # ---___---___---___--- Load and Filter Tree Sequence ---___---___---___--- #

    tree_file = os.path.join(config["output_folder"], "simulation.trees")
    tree_sequence = tskit.load(tree_file)

    kept_individuals = [ind.id for ind in tree_sequence.individuals() if 
                        (ind.flags & pyslim.INDIVIDUAL_REMEMBERED) or 
                        (ind.flags & pyslim.INDIVIDUAL_RETAINED)]
    kept_nodes = [node for ind in kept_individuals for node in tree_sequence.individual(ind).nodes]

    filtered_ts = tree_sequence.simplify(kept_nodes, keep_input_roots=True)

    # ---___---___---___--- Recapitation ---___---___---___--- #

    demography = msprime.Demography()
    demography.add_population(name="p1", initial_size=config["recap_Ne"])
    recap_ts = pyslim.recapitate(filtered_ts, recombination_rate=1e-8, demography=demography)

    # ---___---___---___--- Add Mutations ---___---___---___--- #

    mut_model = msprime.SMM(lo=5, hi=50)
    mut_ts = msprime.sim_mutations(recap_ts, rate=config["mutation_rate"], model=mut_model, random_seed=config["seed"])

    # ---___---___---___--- Format Data for NeEstimator (GENEPOP) ---___---___---___--- #

    def copy_number_matrix(ts):
        C = np.zeros((ts.num_sites, ts.num_samples), dtype=int)
        for var in ts.variants():
            alleles = np.array([int(allele) if allele is not None else 0 for allele in var.alleles])
            C[var.site.id] = alleles[var.genotypes]
        return C

    copy_numbers = copy_number_matrix(mut_ts).T  # shape: (num_individuals, num_loci)
    num_individuals, num_loci = copy_numbers.shape
    output_gen_file = os.path.join(sim_folder, "simulation_data.gen")

    # Separate the different samples and mark "pop" between them
    sample_sizes_Ne = config["sample_sizes_Ne"]
    sample1_size_Ne = sample_sizes_Ne[0]
    sample2_size_Ne = sample_sizes_Ne[1]

    with open(output_gen_file, "w") as f:
        # Title
        f.write("Simulated GENEPOP Data\n")
        f.write(" ".join([f"Locus_{i+1}" for i in range(num_loci)]) + "\n")

        # First sample
        f.write("pop\n")
        for i in range(0, sample1_size_Ne * 2, 2):
            genotype_line = ", " + " ".join(
                f"{int(copy_numbers[i, j]):03}{int(copy_numbers[i+1, j]):03}"
                for j in range(num_loci)
            )
            f.write(f"Indiv_{(i//2)+1}{genotype_line}\n")

        # Second sample
        f.write("pop\n")
        offset = sample1_size_Ne * 2
        for i in range(offset, offset + sample2_size_Ne * 2, 2):
            genotype_line = ", " + " ".join(
                f"{int(copy_numbers[i, j]):03}{int(copy_numbers[i+1, j]):03}"
                for j in range(num_loci)
            )
            f.write(f"Indiv_{(i//2)+1}{genotype_line}\n")

    # ---___---___---___--- Run NeEstimator ---___---___---___--- #

    ne2_exe = os.path.join(SCRIPT_DIR, "Ne2.exe")

    # --- Copy info file template to the simulation file ---
    os.makedirs(os.path.join(SCRIPT_DIR, "templates"), exist_ok=True)
    template_info = os.path.join(SCRIPT_DIR, "templates", "info")

    shutil.copy(template_info, os.path.join(sim_folder, "info"))

    try:
        subprocess.run([ne2_exe, "i:info"], cwd=sim_folder, check=True)
        # Moving the data file created to the right directory
        simulations_dir = os.path.join(SCRIPT_DIR, "simulations")
        origin_file = os.path.join(simulations_dir, "simulation_dataNe.txt")
        destination_file = os.path.join(sim_folder, "simulation_dataNe.txt")

        if os.path.exists(origin_file):
            try:
                shutil.move(origin_file, destination_file)
            except Exception as e:
                pass
        else:
            pass

    except subprocess.CalledProcessError:
        pass

    # ---___---___---___--- Organize the data into a CSV file ---___---___---___--- #

    def extract_ne_stats(txt_path):
        """Extrait les principales stats du fichier de NeEstimator"""

        with open(txt_path, "r") as f:
            content = f.read()

        # ---- LINKAGE / HETERO / COAN ---- (Population 2 uniquement)
        ne_values = re.findall(r"Population\s+2.*?Estimated Ne\^ =.*?([\d.]+)", content, re.DOTALL)
        neb_values = re.findall(r"Population\s+2.*?Estimated Neb\^  =.*?([\d.]+)", content, re.DOTALL)
        neb_coan = re.findall(r"Population\s+2.*?MOLECULAR COANCESTRY METHOD.*?Estimated Neb\^ =\s+([\d.]+|Infinite)", content, re.DOTALL)

        # ---- MÉTHODES TEMPORELLES ----
        # Pollak
        pollak_ne = re.findall(r"\(Pollak\).*?\* Ne =\s+([-]?\d+\.\d+)", content, re.DOTALL)
        # Nei/Tajima
        nei_ne = re.findall(r"\(Nei/Tajima\).*?\* Ne =\s+([-]?\d+\.\d+)", content, re.DOTALL)
        # Jorde/Ryman
        jorde_ne = re.findall(r"\(Jorde/Ryman\).*?\* Ne =\s+([-]?\d+\.\d+)", content, re.DOTALL)

        # Formatage final
        def parse_value(v):
            try:
                f = float(v)
                return f if f > 0 else None
            except:
                return None

        return {
            "Ne_LD": parse_value(ne_values[0]) if len(ne_values) > 0 else None,
            "Neb_HE": round(sum(map(parse_value, neb_values[:3])) / 3, 2) if neb_values else None,
            "Neb_Coan": parse_value(neb_coan[0]) if neb_coan else None,
            "Ne_Pollak": parse_value(pollak_ne[0]) if pollak_ne else None,
            "Ne_Nei": parse_value(nei_ne[0]) if nei_ne else None,
            "Ne_Jorde": parse_value(jorde_ne[0]) if jorde_ne else None, 
        }

    def read_config(path):
        """Lit le fichier slim_config.txt en dict, avec traitement spécial des sample_sizes"""
        config_dict = {}
        with open(path, "r") as f:
            for line in f:
                if "=" in line:
                    key, value = line.strip().split("=", 1)

                    # Cas des tailles d’échantillons structurées
                    if key in ["sample_sizes_Ne", "sample_sizes_CMR"]:
                        try:
                            s1, s2 = map(int, value.split(","))
                            suffix = key.replace("sample_sizes_", "")
                            config_dict[f"sample1_size_{suffix}"] = s1
                            config_dict[f"sample2_size_{suffix}"] = s2
                        except ValueError:
                            pass
                    else:
                        config_dict[key] = value
        return config_dict

    # Fusion et écriture finale
    summary_path = os.path.join(all_simulations, "summary_table.csv")
    config_dict = read_config(slim_config_file)

    # Ajouter les stats de NeEstimator
    ne_data_path = os.path.join(sim_folder, "simulation_dataNe.txt")
    if os.path.exists(ne_data_path):
        ne_stats = extract_ne_stats(ne_data_path)
        config_dict.update(ne_stats)
    else:
        pass

    # Ajouter l'identifiant de simulation
    config_dict["simulation_id"] = sim_id

    # Supprimer les chemins et infos internes
    UNWANTED_KEYS = ["output_folder", "log_file", "timestamp", "seed"]
    for key in UNWANTED_KEYS:
        config_dict.pop(key, None)

    # Enregistrement dans le tableau global
    df_row = pd.DataFrame([config_dict])

    if os.path.exists(summary_path):
        df_existing = pd.read_csv(summary_path)
        df_combined = pd.concat([df_existing, df_row], ignore_index=True)
        df_combined.to_csv(summary_path, index=False)
    else:
        df_row.to_csv(summary_path, index=False)



