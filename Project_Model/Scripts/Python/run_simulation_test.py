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
<<<<<<< Updated upstream
        "pop_size" : 200,
        "num_loci" : 19,
        "num_generations" : 9,
        "sample_sizes_Ne" : [50, 50],
        "sample_sizes_CMR" : [50, 50],
        "mutation_rate" : 5*10e-4,
        "recap_Ne" : 200,
=======
        "pop_size" : 100,
        "num_loci" : 20,
        "num_generations" : 40,
        "sample_sizes_Ne" : [50,50],
        "sample_sizes_CMR" : [100,100],
        "low_repeats" : 1,
        "high_repeats" : 200,
        "mutation_rate" : 1e-2,
        "recap_Ne" : 100,
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
    mut_model = msprime.SMM(lo=5, hi=50)
=======
    lo_repeat = config["low_repeats"]
    hi_repeat = config["high_repeats"]
    root_dist = [0.0] * (hi_repeat - lo_repeat + 1)
    root_dist[100 - hi_repeat - 1] = 1.0 
    mut_model = msprime.SMM(lo=lo_repeat, hi=hi_repeat, root_distribution = root_dist)
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
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
=======
        for pop in [1, 2]:
            # LINKAGE
            ne_values_raw = re.findall(rf"Population\s+{pop}.*?Estimated Ne\^ =\s*(\S+)", content, re.DOTALL)
            ld_ne_value = float(ne_values_raw[0]) if ne_values_raw and ne_values_raw[0] not in ("Infinite", "None") else None

            # HETEROZYGOTE EXCESS (par bloc population)
            he_block = re.search(rf"Population\s+{pop}.*?HETEROZYGOTE EXCESS METHOD.*?Estimated Neb\^  =\s+(\S+)", content, re.DOTALL)
            if he_block:
                he_raw = he_block.group(1)
                he_val = float(he_raw) if "Inf" not in he_raw else None
            else:
                he_val = None

            # COANCESTRY
            coan_block = re.search(rf"Population\s+{pop}.*?MOLECULAR COANCESTRY METHOD.*?Estimated Neb\^ =\s+(\S+)", content, re.DOTALL)
            if coan_block:
                coan_val = float(coan_block.group(1)) if "Inf" not in coan_block.group(1) else None
            else:
                coan_val = None

            # Temporal methods (Two-sample methods)
            if pop == 2:  # Only with two samples
                pollak_ne = re.findall(r"\(Pollak\).*?\* Ne =\s+([-]?\d+\.\d+)", content, re.DOTALL)
                nei_ne = re.findall(r"\(Nei/Tajima\).*?\* Ne =\s+([-]?\d+\.\d+)", content, re.DOTALL)
                jorde_ne = re.findall(r"\(Jorde/Ryman\).*?\* Ne =\s+([-]?\d+\.\d+)", content, re.DOTALL)

            def parse_value(v):
                try:
                    f = float(v)
                    return f if f > 0 else None
                except:
                    return None

            # One-sample methods
            results.update({
                f"LD_Ne_0.05_Pop{pop}": ld_ne_value,
                f"HE_Neb_mean_Pop{pop}": he_val,
                f"Coan_Neb_n_Pop{pop}": coan_val
            })

        # Temporal methods
        results.update({
>>>>>>> Stashed changes
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
    else:
        pass

    # Add the simulation ID
    config_dict["simulation_id"] = sim_id

# Mean allele number for each sample
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


    summary_txt_path = os.path.join(sim_folder, "summary.txt")
    with open(summary_txt_path, "w") as f:
        def write_section(header, keys, file_handle):
            file_handle.write(f"\n[{header}]\n")
            for key in keys:
                if key in config_dict:
                    file_handle.write(f"{key:<24} = {config_dict[key]}\n")

    with open(summary_txt_path, "w") as f:
        write_section("Simulation Info", ["simulation_id", "timestamp", "seed", "output_folder"], f)
        write_section("Model Parameters", ["pop_size", "num_loci", "num_generations", "low_repeats", "high_repeats", "mutation_rate", "recap_Ne"], f)
        write_section("Sampling Design", ["sample1_size_Ne", "sample2_size_Ne", "sample1_size_CMR", "sample2_size_CMR"], f)
        write_section("Population Census", ["census_N"], f)
        write_section("Ne Estimates - One Sample", [
                "LD_Ne_0.05_Pop1", "HE_Neb_mean_Pop1", "Coan_Neb_n_Pop1",
                "LD_Ne_0.05_Pop2", "HE_Neb_mean_Pop2", "Coan_Neb_n_Pop2"
            ], f)
        write_section("Ne Estimates - Temporal", ["Ne_Pollak", "Ne_Nei", "Ne_Jorde"], f)
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


    # Delete unsollicited parts of the config file
    UNWANTED_KEYS = ["output_folder", "log_file", "timestamp", "seed"]
    for key in UNWANTED_KEYS:
        config_dict.pop(key, None)

    # Save the summary file
    df_row = pd.DataFrame([config_dict])

    if os.path.exists(summary_path):
        df_existing = pd.read_csv(summary_path)
        df_combined = pd.concat([df_existing, df_row], ignore_index=True)
        df_combined.to_csv(summary_path, index=False)
    else:
        df_row.to_csv(summary_path, index=False)

        import os

    # Liste des noms de fichiers à supprimer
    files_to_remove = [
        "info",
        "option",
        "simulation_dataBur.txt",
        "simulation_dataNe.txt",
        "simulation_dataNexHt.txt",
        "simulation_dataNexCn.txt",
        "simulation_dataNexLD.txt",
        "simulation_dataNexTp.txt",
        "simulation_dataLoc.txt",
        "simulation_data.gen",
        "simulation.trees",
        "slim_config.txt",
        "slim.log"
    ]

    # Suppression
    for filename in files_to_remove:
        filepath = os.path.join(sim_folder, filename)
        if os.path.exists(filepath):
            try:
                os.remove(filepath)
            except Exception as e:
                pass




