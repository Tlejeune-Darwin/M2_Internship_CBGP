import re
import numpy as np #type: ignore

def clean_value(v):
    try:
        if isinstance(v, str) and v.strip().lower() in ("infinite", "inf", "none", ""):
            return None
        return float(v)
    except:
        return None

def extract_values(label, text):
    match = re.search(rf"{label}\s*=\s*(.+)", text)
    if match:
        return [clean_value(v) for v in re.split(r"\s{2,}", match.group(1).strip())]
    return [None] * 4

def parse_ld(content, pop):
    results = {}
    block = re.search(rf"Population\s+{pop}.*?LINKAGE DESEQUILIBRIUM METHOD.*?HETEROZYGOTE EXCESS METHOD", content, re.DOTALL)
    if block:
        text = block.group(0)
        ne_vals = extract_values("Estimated Ne\\^", text)
        r2_vals = extract_values("OverAll r\\^2", text)
        for i, th in enumerate(["0.05", "0.02", "0.01", "0.00"]):
            results[f"LD_Ne_{th}_Pop{pop}"] = ne_vals[i]
            results[f"r2_overall_{th}_Pop{pop}"] = r2_vals[i]
        results[f"LD_Ne_Pop{pop}"] = ne_vals
        results[f"LD_r2_Pop{pop}"] = r2_vals
    return results

def parse_he(content, pop):
    results = {}
    block = re.search(rf"Population\s+{pop}.*?HETEROZYGOTE EXCESS METHOD.*?MOLECULAR COANCESTRY METHOD", content, re.DOTALL)
    if block:
        text = block.group(0)
        he_vals = extract_values("Estimated Neb\\^", text)
        d_vals = extract_values("Weighted Mean D", text)
        for i, th in enumerate(["0.05", "0.02", "0.01", "0.00"]):
            results[f"HE_Neb_mean_{th}_Pop{pop}"] = he_vals[i]
            results[f"HE_weighted_D_mean_{th}_Pop{pop}"] = d_vals[i]
        results[f"HE_Neb_mean_Pop{pop}"] = he_vals
        results[f"HE_Neb_weighted_D_mean_Pop{pop}"] = d_vals
    return results

def parse_coancestry(content, pop):
    results = {}
    block = re.search(rf"Population\s+{pop}.*?MOLECULAR COANCESTRY METHOD.*?Estimated Neb\\^ =\s+(\S+)", content, re.DOTALL)
    f1_block = re.search(rf"Population\s+{pop}.*?MOLECULAR COANCESTRY METHOD.*?OverAll f1\\^.*?=\s+(-?\d+\.\d+)", content, re.DOTALL)
    results[f"Coan_Neb_n_Pop{pop}"] = clean_value(block.group(1)) if block else None
    results[f"Coan_f1_Pop{pop}"] = clean_value(f1_block.group(1)) if f1_block else None
    return results

def parse_pollak(content):
    results = {"P_Ne": [None]*4, "P_Fk": [None]*4, "P_F'": [None]*4}
    block = re.search(r"\(Pollak\)(.*?)(?=\(Nei/Tajima\)|\Z)", content, re.DOTALL)
    if block:
        text = block.group(1)
        results["P_Ne"] = extract_values(r"\* Ne", text)
        results["P_Fk"] = extract_values("Fk", text)
        results["P_F'"] = extract_values("F'", text)
    return results

def parse_nei(content):
    results = {"N_Ne": [None]*4, "N_Fc": [None]*4, "N_F'": [None]*4}
    block = re.search(r"\(Nei/Tajima\)(.*?)(?=\(Jorde/Ryman\)|\Z)", content, re.DOTALL)
    if block:
        text = block.group(1)
        results["N_Ne"] = extract_values(r"\* Ne", text)
        results["N_Fc"] = extract_values("Fc", text)
        results["N_F'"] = extract_values("F'", text)
    return results

def parse_jorde(content):
    results = {"J_Ne": [None]*4, "J_Fs": [None]*4, "J_F'": [None]*4}
    block = re.search(r"\(Jorde/Ryman\)(.*?)(?=Ending time:|\Z)", content, re.DOTALL)
    if block:
        text = block.group(1)
        results["J_Ne"] = extract_values(r"\* Ne", text)
        results["J_Fs"] = extract_values("Fs", text)
        results["J_F'"] = extract_values("F'", text)
    return results

def extract_ne_stats(txt_path):
    with open(txt_path, "r") as f:
        content = f.read()

    results = {}
    for pop in [1, 2]:
        results.update(parse_ld(content, pop))
        results.update(parse_he(content, pop))
        results.update(parse_coancestry(content, pop))

    results.update(parse_pollak(content))
    results.update(parse_nei(content))
    results.update(parse_jorde(content))

    return results
