def run_simulation_linux(base_dir="simulations", pop_size=None, num_loci=None, sim_prefix="sim"):

    # ---___---___---___--- 1. Imports ---___---___---___--- #
    
    # Packages needed to run python script
    import os
    import sys
    import subprocess
    import tskit                            # type: ignore
    import pyslim                           # type: ignore
    import msprime                          # type: ignore
    import numpy as np                      # type: ignore
    import random
    import warnings
    import shutil
    import pandas as pd                     # type: ignore
    import re
    from datetime import datetime
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    if SCRIPT_DIR not in sys.path:
        sys.path.append(SCRIPT_DIR)
    from labels import better_names
    inverse_better_names = {v: k for k, v in better_names.items()}
    
    # ---___---___---___--- 2. Initialization and Paths ---___---___---___--- #

    ### 2.1. Create a global config file for later study purpose ###
    def get_global_config(simulations_dir):
        """Vérifie si le fichier config_global.txt existe, le crée sinon, puis le lit"""
        global_config_path = os.path.join(simulations_dir, "config_global.txt")
        
        # Create the file if not done already
        if not os.path.exists(global_config_path):
            with open(global_config_path, "w") as f:
                f.write("# Configuration générale des simulations\n")
                f.write(f"{better_names['num_loci']} = 20\n")
                f.write(f"{better_names['low_repeats']} = 1\n")
                f.write(f"{better_names['high_repeats']} = 200\n")
                f.write(f"{better_names['mutation_rate']} = 0.00001, 0.001\n")
                f.write(f"{better_names['sample1_generation']} = 30\n")
                f.write(f"{better_names['sample2_generation']} = 10\n")
                f.write(f"{better_names['sample_sizes_Ne']} = 50,50\n")
                f.write(f"{better_names['sample_sizes_CMR']} = 100\n")
                f.write(f"{better_names['pop_size_logrange']} = 50,10000\n")

        # Reading the file
        config = {}
        with open(global_config_path, "r") as f:
            for line in f:
                if "=" in line and not line.strip().startswith("#"):
                    raw_key, value = line.strip().split("=", 1)
                    key = inverse_better_names.get(raw_key.strip(), raw_key.strip())
                    config[key] = value.strip()
        return config

    ### 2.2. Directory script ###
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    
    ### 2.3. Choose depending on the desktop name ###
    def get_desktop_path():
        desktop_fr = os.path.join(os.path.expanduser("~"), "Bureau")
        desktop_en = os.path.join(os.path.expanduser("~"), "Desktop")
        return desktop_fr if os.path.isdir(desktop_fr) else (desktop_en if os.path.isdir(desktop_en) else os.path.expanduser("~"))

    ### 2.4. Simulations directory placed on the desktop ###
    all_simulations = base_dir
    os.makedirs(all_simulations, exist_ok=True)
    
    ### 2.5. Create a folder for each simulation ###
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    existing = [d for d in os.listdir(all_simulations) if d.startswith(sim_prefix)]
    numbers = [int(re.search(rf"{sim_prefix}_(\d+)", name).group(1)) for name in existing if re.search(rf"{sim_prefix}_(\d+)", name)]
    next_sim_num = max(numbers, default=0) + 1

    # Format : sim_0000001 up to 1 milllion
    sim_id = f"{sim_prefix}_{next_sim_num:07d}"
    sim_folder = os.path.join(all_simulations, sim_id)
    os.makedirs(sim_folder, exist_ok=True)


    # ---___---___---___--- 3. Config File Generation ---___---___---___--- #

    ### 3.1. Create the simulation parameter dictionary "config_file" ###
    global_config = get_global_config(all_simulations)
    if pop_size is None:
        low, high = map(float, global_config["pop_size_logrange"].split(","))
        pop_size = int(np.exp(np.random.uniform(np.log(low), np.log(high))))
        log_mu_low = np.log(1e-5)
        log_mu_high = np.log(1e-3)
        mutation_rate = float(np.exp(np.random.uniform(log_mu_low, log_mu_high)))

    config = {
        "simulation_id" : sim_id,                                                               # Name of this specific simulation
        "pop_size" : pop_size,                                                                  # Size of the population
        "num_loci" : int(global_config["num_loci"]),                                            # Number of loci used in SLiM   
        "sample1_generation" : int(global_config["sample1_generation"]),                        # Number n of generations before the first genetic sample is taken
        "sample2_generation" : int(global_config["sample2_generation"]),                        # Number n of generations between the two genetic samples
        "sample_sizes_Ne" : list(map(int, global_config["sample_sizes_Ne"].split(","))),        # Size in individuals of the genetic samples   
        "sample_sizes_CMR" : int(global_config["sample_sizes_CMR"]),                            # Size in individuals of the demographic samples
        "low_repeats" : int(global_config["low_repeats"]),                                      # Lowest number of repeats for simulation of mutation processes
        "high_repeats" : int(global_config["high_repeats"]),                                    # Highest number of repeats
        "mutation_rate" : mutation_rate,                                                        # Mutation rate used during mutation simulation
        "recap_Ne" : pop_size,                                                                  # Effective size attributed for the recapitation
        "output_folder" : sim_folder,                                                           # All simulations end up in the main sim_folder
        "timestamp" : timestamp,                                                                # Timestamp placed in the name of each simulation
        "seed" : random.randint(1, 10**6)                                                       # For analysis purpose
    }

    ### 3.2. Write the "slim_config.txt" file (for SLiM input) ###
    slim_config_file = os.path.join(sim_folder, "slim_config.txt")
    with open(slim_config_file, "w") as f:
        for key, value in config.items():
            if isinstance(value, list):
                value = ",".join(map(str, value))
            f.write(f"{key}={value}\n")

    ### 3.3. Create the "info" file (for NeEstimator input) ###
    info_path = os.path.join(sim_folder, "info")
    sample1_generation = int(config["sample1_generation"])
    sample2_generation = int(config["sample2_generation"])

    with open(info_path, "w") as f:
        f.write("15  0\n")                                                                      # Calculation methods : by addition (15 means all methods)
        f.write("\n")                                                                           # The line break is important here
        f.write("simulation_data.gen\n")                                                        # Input directory
        f.write("2\n")                                                                          # Define the format of the file (1 for FSTAT; 2 for GENEPOP)
        f.write(os.path.join(sim_folder, "") + "\n")                                            # Output directory
        f.write("simulation_dataNe.txt\n")                                                      # Output folder name
        f.write("3\n")                                                                          # Number of critical values
        f.write("0.05  0.02  0.01\n")                                                           # Critical values
        f.write("0\n")                                                                          # Random mating (0) or monogamy (1)
        f.write(f"0 {sample1_generation} {sample2_generation+sample1_generation}\n")            # 3 numbers here : first one represents the population size if we are in Plan I, second one represents the generation of the first sample and last one represents the generation of the second sample
        f.write("0\n")                                                                          # 0 must be added here to stop the generations

    ### 3.4. Create the "option" file (for NeEstimator input) ###
    option_path = os.path.join(sim_folder, "option")

    option_lines = [
        "15  0  1  1",                                                                          # First number : sum of methods ; second number : sum of temporal methods (7 = all TP methods); third number : number of critical values added (if 0 only the smallest is shown); Fourth entry : tab delimiter (if > 0)
        "0",                                                                                    # Maximum individuals/pop, if 0 : no limit
        "-1",                                                                                   # -1 to ouptut the allelic frequencies
        "-1  1  0  0",                                                                          # -1 activate Burrows outputs for all pop; second entry shows all critical values
        "0",                                                                                    # Parameter Confidence Interval (Deactivated)
        "0",                                                                                    # Jackknife Confidence Interval (Deactivated)
        "0",                                                                                    # Up to population, or range of population to run (if more than 2), 0 means no restriction
        "0",                                                                                    # All loci are accepted
        "1",                                                                                    # Input 1 to create a file that documents missing data from the input file
        "0"                                                                                     # Line for chromosome/Loci option and file
    ]

    # Alternative method of writing in the file (other than the one used for the "info" file)
    with open(option_path, "w") as f:
        f.write("\n".join(option_lines) + "\n")

    # ---___---___---___--- 4. Run SLiM Simulation ---___---___---___--- #

    ### 4.1. Define paths to the SLiM executable and script ###
    slim_script = os.path.join(SCRIPT_DIR, "..", "SLiM", "Sim_model.slim")
    slim_executable = os.path.join(SCRIPT_DIR, "..", "..", "Bin", "slim")
    slim_config_file = slim_config_file.replace("\\", "/")

    ### 4.2. Build the SLiM command ###
    slim_command = [slim_executable, "-d", f'config_file="{slim_config_file}"', slim_script]

    log_file_path = os.path.join(sim_folder, "slim.log")

    ### 4.3. Run SLiM via subprocess and log output to slim.log ###
    with open(log_file_path, "w") as log_file:
        slim_process = subprocess.run(slim_command, stdout = log_file, stderr = log_file, text = True, cwd = SCRIPT_DIR)

    # Delete the TimeUnitMismatch Warning 
    warnings.simplefilter("ignore", msprime.TimeUnitsMismatchWarning)

    # 4.4. Lire les valeurs de MatchCount_i et Census_N_i
    match_census_path = os.path.join(sim_folder, "census_match_output.txt")
    if os.path.exists(match_census_path):
        with open(match_census_path, "r") as f:
            for line in f:
                if "=" in line:
                    key, val = line.strip().split("=")
                    try:
                        config[key] = int(val)
                    except ValueError:
                        config[key] = val


    # ---___---___---___--- 5. Tree Sequence Processing ---___---___---___--- #

    ### 5.1. Load the ".trees" file ###
    tree_file = os.path.join(config["output_folder"], "simulation.trees")
    tree_sequence = tskit.load(tree_file)

    ### 5.2. Filter for "REMEMBERED" and "RETAINED" individuals ###
    kept_individuals = [ind.id for ind in tree_sequence.individuals() if 
                        (ind.flags & pyslim.INDIVIDUAL_REMEMBERED) or 
                        (ind.flags & pyslim.INDIVIDUAL_RETAINED)]
    kept_nodes = [node for ind in kept_individuals for node in tree_sequence.individual(ind).nodes]

    ### 5.3. Simplify the tree sequence to reduce memory and noise ###
    filtered_ts = tree_sequence.simplify(kept_nodes, keep_input_roots=True)

    # ---___---___---___--- 6. Recapitation ---___---___---___--- #

    ### 6.1. Create the demographic model ###
    demography = msprime.Demography()
    demography.add_population(name="p1", initial_size=config["recap_Ne"])
    
    ### 6.2. Recapitate ###
    recap_ts = pyslim.recapitate(filtered_ts, recombination_rate=1e-8, demography=demography)

    # ---___---___---___--- 7. Simulate Mutations ---___---___---___--- #

    ### 7.1. Define the stepwise mutation model (SMM) ###
    lo_repeat = config["low_repeats"]
    hi_repeat = config["high_repeats"]
    root_dist = [0.0] * (hi_repeat - lo_repeat + 1)
    root_dist[100 - hi_repeat - 1] = 1.0 
    mut_model = msprime.SMM(lo=lo_repeat, hi=hi_repeat, root_distribution = root_dist)

    ### 7.2. Simulate mutation process ###
    mut_ts = msprime.sim_mutations(recap_ts, rate=config["mutation_rate"], model=mut_model, random_seed=config["seed"])

    # ---___---___---___--- 8. Format Data for NeEstimator (GENEPOP) ---___---___---___--- #

    ### 8.1. Generate the copy number matrix ###
    def copy_number_matrix(ts):
        C = np.zeros((ts.num_sites, ts.num_samples), dtype=int)
        for var in ts.variants():
            alleles = np.array([int(allele) if allele is not None else 0 for allele in var.alleles])
            C[var.site.id] = alleles[var.genotypes]
        return C

    copy_numbers = copy_number_matrix(mut_ts).T  # shape: (num_individuals, num_loci)
    num_individuals, num_loci = copy_numbers.shape
    output_gen_file = os.path.join(sim_folder, "simulation_data.gen")

    ### Write the ".gen" file in GENEPOP format ###
    sample_sizes_Ne = config["sample_sizes_Ne"]
    sample2_size_Ne = sample_sizes_Ne[0]
    sample1_size_Ne = sample_sizes_Ne[1]

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

    # ---___---___---___--- 9. Run NeEstimator ---___---___---___--- #

    ### 9.1. Execute the NeEstimator binary (Ne2x)
    ne2_exe = os.path.join(SCRIPT_DIR, "..", "..", "Bin", "Ne2x")

    try:
        subprocess.run([ne2_exe, "i:info", "o:option"], cwd=sim_folder, check=True)
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

    # ---___---___---___--- 10. Extract NeEstimator results ---___---___---___--- #

    def extract_ne_stats(txt_path):
        """Extract the main stats from the NeEstimator analysis"""

        with open(txt_path, "r") as f:
            content = f.read()

            def clean_value(v):
                try:
                    if isinstance(v, str) and v.strip().lower() in ("infinite", "inf", "none", ""):
                        return None
                    return float(v)
                except:
                    return None

        results = {}

        ### 10.1. Parse One-sample estimates : LD, HE, Coancestry ###
        for pop in [1, 2]:
            ne_vals = [None] * 4
            r2_vals = [None] * 4
            he_vals = [None] * 4 
            d_vals = [None] * 4
            coan_val = None
            f1_val = None

            results[f"LD_Ne_Pop{pop}"] = [None] * 4
            results[f"LD_r2_Pop{pop}"] = [None] * 4
            results[f"He_Neb_mean_Pop{pop}"] = [None] * 4
            results[f"He_weighted_D_mean_Pop{pop}"] = [None] * 4
            results[f"Coan_Neb_n_Pop{pop}"] = None
            results[f"Coan_f1_Pop{pop}"] = None

        results.update({
            "P_Ne": [None] * 4,
            "P_Fk": [None] * 4,
            "P_Fprime": [None] * 4,
            "N_Ne": [None] * 4,
            "N_Fc": [None] * 4,
            "N_Fprime": [None] * 4,
            "J_Ne": [None] * 4,
            "J_Fs": [None] * 4,
            "J_Fprime": [None] * 4,
        })
            
        for pop in [1, 2]:

            # LINKAGE DESEQUILIBRIUM
            ld_block = re.search(rf"Population\s+{pop}.*?LINKAGE DISEQUILIBRIUM METHOD.*?HETEROZYGOTE EXCESS METHOD", content, re.DOTALL)
            if ld_block:
                ld_text = ld_block.group(0)

                def extract_four_values(label):
                    match = re.search(rf"{label}\s*=\s*(.+)", ld_text)
                    if match:
                        values = re.split(r"\s{2,}", match.group(1).strip())
                        return [float(v) if v != "Infinite" else None for v in values]
                    return [None] * 4

                ne_vals = extract_four_values("Estimated Ne\\^")
                r2_vals = extract_four_values("OverAll r\\^2")

            # HETEROZYGOTE EXCESS
            # === HETEROZYGOTE EXCESS PARSE ===
            he_block = re.search(rf"Population\s+{pop}.*?HETEROZYGOTE EXCESS METHOD.*?MOLECULAR COANCESTRY METHOD", content, re.DOTALL)
            if he_block:
                he_text = he_block.group(0)

                # Extraction similaire à LD
                def extract_he_values(label):
                    match = re.search(rf"{label}\s*=\s*(.+)", he_text)
                    if match:
                        values = re.split(r"\s{2,}", match.group(1).strip())
                        return [float(v) if v != "Infinite" else None for v in values]
                    return [None] * 4

                he_vals = extract_he_values("Estimated Neb\\^")
                d_vals = extract_he_values("Weighted Mean D")

            # COANCESTRY
            coan_block = re.search(rf"Population\s+{pop}.*?MOLECULAR COANCESTRY METHOD.*?Estimated Neb\^ =\s+(\S+)", content, re.DOTALL)
            if coan_block:
                coan_val = float(coan_block.group(1)) if "Inf" not in coan_block.group(1) else None
            else:
                coan_val = None
            f1_block = re.search(rf"Population\s+{pop}.*?MOLECULAR COANCESTRY METHOD.*?OverAll f1\^.*?=\s+(-?\d+\.\d+)", content, re.DOTALL)
            if f1_block:
                f1_val = float(f1_block.group(1))
            else:
                f1_val = None
            
            ### 10.2. Parse Two-sample estimates : Pollak, Nei, Jorde ###
            # Temporal methods (Two-sample methods)
            if pop == 2:  # Only with two samples
                def extract_temporal_values(label, text):
                    match = re.search(rf"{label}\s*=\s*(.+)", text)
                    if match:
                        values = re.split(r"\s{2,}", match.group(1).strip())
                        return [clean_value(v) for v in values]
                    return [None] * 4

                # Pollak Method
                pollak_vals = [None] * 4
                fk_vals = [None] * 4
                P_fprime_vals = [None] * 4
                pollak_block = re.search(r"\(Pollak\)(.*?)(?=\(Nei/Tajima\)|\Z)", content, re.DOTALL)
                if pollak_block:
                    pollak_text = pollak_block.group(1)
                    pollak_vals = extract_temporal_values(r"\* Ne", pollak_text)

                    fk_match = re.search(r"Fk\s*=\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", pollak_text)
                    if fk_match:
                        fk_vals = [clean_value(fk_match.group(i)) for i in range(1, 5)]

                    fprime_match = re.search(r"F'\s*=\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)", pollak_text)
                    if fprime_match:
                        P_fprime_vals = [clean_value(fprime_match.group(i)) for i in range(1, 5)]

                # Nei/Tajima Method
                nei_vals = [None] * 4
                fc_vals = [None] * 4
                N_fprime_vals = [None] * 4
                nei_block = re.search(r"\(Nei/Tajima\)(.*?)(?=\(Jorde/Ryman\)|\Z)", content, re.DOTALL)
                if nei_block:
                    nei_text = nei_block.group(1)
                    nei_vals = extract_temporal_values(r"\* Ne", nei_text) if nei_text else [None]*4

                    fc_match = re.search(r"Fc\s*=\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", nei_text)
                    if fc_match:
                        fc_vals = [clean_value(fc_match.group(i)) for i in range(1, 5)]

                    N_fprime_match = re.search(r"F'\s*=\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)", nei_text)
                    if N_fprime_match:
                        N_fprime_vals = [clean_value(N_fprime_match.group(i)) for i in range(1, 5)]

                # Jorde/Ryman Method
                jorde_vals = [None] * 4
                fs_vals = [None] * 4
                J_fprime_vals = [None] * 4
                jorde_block = re.search(r"\(Jorde/Ryman\)(.*?)(?=Ending time:|\Z)", content, re.DOTALL)
                if jorde_block:
                    jorde_text = jorde_block.group(1)
                    jorde_vals = extract_temporal_values(r"\* Ne", jorde_text) if jorde_text else [None]*4

                    fs_match = re.search(r"Fs\s*=\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", jorde_text)
                    if fs_match:
                        fs_vals = [clean_value(fs_match.group(i)) for i in range(1, 5)]

                    J_fprime_match = re.search(r"F'\s*=\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)", jorde_text)
                    if J_fprime_match:
                        J_fprime_vals = [clean_value(J_fprime_match.group(i)) for i in range(1, 5)]

            ### 10.3. Use "parse_value" for cleanup and validation ###
            def parse_value(v):
                try:
                    f = float(v)
                    return f if f > 0 else None
                except:
                    return None

            # One-sample methods
            results.update({
                f"LD_Ne_Pop{pop}": ne_vals,
                f"LD_r2_Pop{pop}": r2_vals,
                f"HE_Neb_mean_Pop{pop}": he_vals,
                f"HE_weighted_D_mean_Pop{pop}" : d_vals,
                f"Coan_Neb_n_Pop{pop}": coan_val,
                f"Coan_f1_Pop{pop}" : f1_val
            })

        # Temporal methods
        results.update({
            "P_Ne": pollak_vals,
            "P_Fk" : fk_vals,
            "P_Fprime" : P_fprime_vals,
            "N_Ne": nei_vals,
            "N_Fc": fc_vals,
            "N_Fprime": N_fprime_vals,
            "J_Ne": jorde_vals,
            "J_Fs": fs_vals,
            "J_Fprime": J_fprime_vals
        })

        return results

    ### 10.4. Load config values with "read_config()" ###
    def read_config(path):
        """Read the config file"""
        config_dict = {}
        with open(path, "r") as f:
            for line in f:
                if "=" in line:
                    key, value = line.strip().split("=", 1)

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

    # Merge the new data
    summary_path = os.path.join(all_simulations, "summary_table.csv")
    config_dict = read_config(slim_config_file)
    ne_data_path = os.path.join(sim_folder, "simulation_dataNe.txt")
    if os.path.exists(ne_data_path):
        ne_stats = extract_ne_stats(ne_data_path)
        config_dict.update(ne_stats)
    thresholds = ["0.050", "0.020", "0.010", "0.000"]
    for pop in [1,2]:
            for label in ["LD_Ne", "LD_r2", "HE_Neb_mean", "HE_weighted_D_mean"]:
                key = f"{label}_Pop{pop}"
                if key in config_dict:
                    values = config_dict.get(key,[])
                    if isinstance(values, list):
                        for i, th in enumerate(thresholds):
                            if i < len(values):
                                config_dict[f"{label}_{th}_Pop{pop}"] = config_dict[key][i]
    
    for label in ["P_Ne", "P_Fk", "P_Fprime", "N_Ne", "N_Fc", "N_Fprime", "J_Ne", "J_Fs", "J_Fprime"]:
        if label in config_dict and isinstance(config_dict[label], list):
            for i, th in enumerate(thresholds):
                if i < len(config_dict[label]):
                    config_dict[f"{label}_{th}"] = config_dict[label][i]
                else:
                    config_dict[f"{label}_{th}"] = None
    else:
        pass

    # Add the simulation ID
    config_dict["simulation_id"] = sim_id

    # ---___---___---___--- 11. Genetic diversity summaries ---___---___---___--- #

    ### 11.1. Read the "simulation_dataLoc.txt" file ###
    loc_file = os.path.join(sim_folder, "simulation_dataLoc.txt")

    per_pop_summaries = {}
    allele_counts = {}
    allele_variances = {}
    het_obs = {}
    het_exp = {}

    if os.path.exists(loc_file):
        with open(loc_file, "r") as f:
            current_pop = None
            locus_counter = 0
            allele_sizes = []

            for line in f:
                match = re.search(r"POPULATION\s+(\d+)", line)
                if match:
                    current_pop = int(match.group(1))
                    per_pop_summaries[current_pop] = []
                    allele_counts[current_pop] = []
                    allele_variances[current_pop] = []
                    het_obs[current_pop] = []
                    het_exp[current_pop] = []
                    locus_counter = 0
                    continue
                
                ### 11.2. Compute per-locus statistics ###
                if re.match(r"\s*\d+:Locus_", line) and current_pop is not None:
                    locus_counter += 1

                elif line.strip().startswith("Alleles:") and current_pop is not None:
                    allele_sizes = list(map(int, line.strip().split("Alleles:")[1].strip().split()))
                    allele_counts[current_pop].append(len(allele_sizes))

                elif line.strip().startswith("Frequencies:") and current_pop is not None:
                    freqs = list(map(float, line.strip().split("Frequencies:")[1].strip().split()))
                    if len(freqs) == len(allele_sizes):
                        mean = sum(a * f for a, f in zip(allele_sizes, freqs))
                        var = sum(f * (a - mean) ** 2 for a, f in zip(allele_sizes, freqs))
                        allele_variances[current_pop].append(var)
                        per_pop_summaries[current_pop].append({
                            "Locus": locus_counter,
                            "Num_Alleles": len(allele_sizes),
                            "Alleles": allele_sizes.copy(),
                            "Frequencies": freqs.copy(),
                            "Variance": round(var, 4)
                        })

                het_match = re.search(r"Overall Het.*?Obs\. \(O\)\s*=\s*([\d.]+),\s*Exp\.?\(E\)\s*=\s*([\d.]+)", line)
                if het_match and current_pop is not None:
                    obs = float(het_match.group(1))
                    exp = float(het_match.group(2))
                    het_obs[current_pop].append(obs)
                    het_exp[current_pop].append(exp)
                    if per_pop_summaries[current_pop]:
                        per_pop_summaries[current_pop][-1]["Obs_Het"] = obs
                        per_pop_summaries[current_pop][-1]["Exp_Het"] = exp

    ### 11.3. Compute per-population averages ###
    for pop_id in per_pop_summaries:
        if het_obs[pop_id]:
            config_dict[f"mean_obs_het_pop{pop_id}"] = round(np.mean(het_obs[pop_id]), 4)
        if het_exp[pop_id]:
            config_dict[f"mean_exp_het_pop{pop_id}"] = round(np.mean(het_exp[pop_id]), 4)
        if allele_counts[pop_id]:
            config_dict[f"mean_alleles_pop{pop_id}"] = round(np.mean(allele_counts[pop_id]), 2)
            config_dict[f"var_alleles_pop{pop_id}"] = round(np.var(allele_counts[pop_id], ddof=1), 2)
            config_dict[f"sum_alleles_pop{pop_id}"] = sum(allele_counts[pop_id])
        if allele_variances[pop_id]:
            config_dict[f"var_allele_size_pop{pop_id}"] = round(np.mean(allele_variances[pop_id]), 4)

    # ---___---___---___--- 12. Write summary ".txt" file ---___---___---___--- #


    ### 12.1. Write overall sections to a human-readable summary ###
    summary_txt_path = os.path.join(sim_folder, "summary.txt")
    with open(summary_txt_path, "w") as f:
        def write_section(header, keys, file_handle):
            file_handle.write(f"\n[{header}]\n")
            filtered_keys = [k for k in keys if k in config_dict]
            if not filtered_keys :
                return
            max_key_len = max(len(k) for k in filtered_keys)
            for key in filtered_keys:
                file_handle.write(f"{key:<{max_key_len}} = {config_dict[key]}\n")

    with open(summary_txt_path, "w") as f:
        write_section("Simulation Info", ["simulation_id", "timestamp", "seed", "output_folder"], f)
        write_section("Model Parameters", ["pop_size", "num_loci", "sample1_generation", "sample2_generation", "low_repeats", "high_repeats", "mutation_rate", "recap_Ne"], f)
        write_section("Sampling Design", ["sample1_size_Ne", "sample2_size_Ne", "sample1_size_CMR", "sample2_size_CMR"], f)
        write_section("Capture-Marquage-Recapture", ["census_N", "MatchCount"], f)
        census_value = config_dict.get("census_N")
        try:
            if census_value is not None and float(census_value) == 0:
                f.write("# There were no match between sampling of individuals during CMR.\n")
        except ValueError:
            pass 

        write_section("Ne Estimates - One Sample - Decreasing critical values [0.050, 0.020, 0.010, 0+]", [], f)
        write_section("Linkage Desequilibrium", ["LD_Ne_Pop1", "LD_r2_Pop1",
                                                 "LD_Ne_Pop2", "LD_r2_Pop2"], f)
        write_section("Heterozygote excess", ["HE_Neb_mean_Pop1", "HE_weighted_D_mean_Pop1",
                                              "HE_Neb_mean_Pop2", "HE_weighted_D_mean_Pop2"], f)
        write_section("Molecular Coancestry", ["Coan_Neb_n_Pop1", "Coan_f1_Pop1",
                                               "Coan_Neb_n_Pop2", "Coan_f1_Pop2"], f)

        write_section("Ne Estimates - Temporal - Decreasing critical values (0.050, 0.020, 0.010, 0+)", [], f)
        write_section("Pollak", ["P_Ne", "P_Fk", "P_Fprime"], f)
        write_section("Nei/Tajima", ["N_Ne", "N_Fc", "N_Fprime"], f)
        write_section("Jorde/Ryman", ["J_Ne", "J_Fs", "J_Fprime"], f)
        
        write_section("Genetic Diversity - Heterozygosity", [
                "mean_exp_het_pop1", "mean_obs_het_pop1",
                "mean_exp_het_pop2", "mean_obs_het_pop2"
            ], f)
        write_section("Genetic Diversity - Allele Counts", [
                "mean_alleles_pop1", "var_alleles_pop1", "sum_alleles_pop1",
                "mean_alleles_pop2", "var_alleles_pop2", "sum_alleles_pop2"
            ], f)
        write_section("Genetic Diversity - Allelic Size Variance", [
                "var_allele_size_pop1", "var_allele_size_pop2"
            ], f)
        
        ### 12.2. Write per-locus summary tables ###
        for pop_id, allelic_summary in per_pop_summaries.items():
            if not allelic_summary:
                continue

            f.write(f"\n[Per-locus summary - Population {pop_id}]\n")

            header_labels = ["Locus", "Num_Alleles", "Alleles", "Frequencies", "Variance", "Obs_Het", "Exp_Het"]

            rows = []
            for idx, entry in enumerate(allelic_summary, start=1):
                row = {
                    "Locus": str(idx),
                    "Num_Alleles": str(entry["Num_Alleles"]),
                    "Alleles": "[" + ", ".join(map(str, entry["Alleles"])) + "]",
                    "Frequencies": "[" + ", ".join(f"{x:.4f}" for x in entry["Frequencies"]) + "]",
                    "Variance": f"{entry['Variance']:.4f}",
                    "Obs_Het": f"{entry.get('Obs_Het', 'NA'):.4f}" if "Obs_Het" in entry else "NA",
                    "Exp_Het": f"{entry.get('Exp_Het', 'NA'):.4f}" if "Exp_Het" in entry else "NA"
                }
                rows.append(row)

            column_widths = {
                key: max(len(key), max(len(row[key]) for row in rows)) for key in header_labels
            }

            header_line = " | ".join(f"{key:<{column_widths[key]}}" for key in header_labels)
            separator_line = "-+-".join("-" * column_widths[key] for key in header_labels)

            f.write(header_line + "\n")
            f.write(separator_line + "\n")
            for row in rows:
                f.write(" | ".join(f"{row[key]:<{column_widths[key]}}" for key in header_labels) + "\n")

    # ---___---___---___--- 13. Append to global CSV summary ---___---___---___--- #

    ### 13.1. Remove unnecessary keys before writing ###
    UNWANTED_KEYS = ["output_folder", "log_file", "timestamp", "seed", "HE_Neb_mean_Pop1", "HE_Neb_mean_Pop2", "HE_weighted_mean_Pop1", "HE_weighted_mean_Pop2"]
    for key in UNWANTED_KEYS:
        config_dict.pop(key, None)
    keys_to_remove = []
    for key, value in config_dict.items():
        if isinstance(value, list) and all(v is None for v in value):
            keys_to_remove.append(key)
    for key in keys_to_remove:
        config_dict.pop(key)

    # Delete the key lists for the csv port (for lisibility)
    clefs_listes = [
        f"{label}_Pop{pop}"
        for label in ["LD_Ne", "LD_r2", "HE_Neb_mean", "HE_weighted_D_mean"]
        for pop in [1, 2]
    ] + ["P_Ne", "P_Fk", "P_Fprime", "N_Ne", "N_Fc", "N_Fprime", "J_Ne", "J_Fs", "J_Fprime"]

    for key in clefs_listes:
        config_dict.pop(key, None)

    ### 13.2. Append the current simulation to "summary_table.csv" ###
    df_row = pd.DataFrame([config_dict])

    if os.path.exists(summary_path):
        df_existing = pd.read_csv(summary_path)
        df_combined = pd.concat([df_existing, df_row], ignore_index=True)
        df_combined.to_csv(summary_path, index=False)
    else:
        df_row.to_csv(summary_path, index=False)

    # ---___---___---___--- 14. Cleanup temporary files ---___---___---___--- #

    ### 14.1. List of files to remove to save space ###
    files_to_remove = [
        "info",
        "option",
        "simulation_dataBur.txt",
        #"simulation_dataNe.txt",
        "simulation_dataNexHt.txt",
        "simulation_dataNexCn.txt",
        "simulation_dataNexLD.txt",
        "simulation_dataNexTp.txt",
        "simulation_dataLoc.txt",
        "simulation_data.gen",
        "simulation.trees",
        "slim_config.txt",
        #"slim.log"
    ]

    ### 14.2. Delete each file if it exists ###
    for filename in files_to_remove:
        filepath = os.path.join(sim_folder, filename)
        if os.path.exists(filepath):
            try:
                os.remove(filepath)
            except Exception as e:
                pass






