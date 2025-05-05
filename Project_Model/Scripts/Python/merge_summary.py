import os, glob
import pandas as pd #type: ignore

base_results_dir = os.path.expanduser("~/results")
batch_dirs = sorted([d for d in os.listdir(base_results_dir) if d.startswith("batch_")])

all_data = []
for batch in batch_dirs:
    batch_path = os.path.join(base_results_dir, batch)
    summary_files = glob.glob(os.path.join(batch_path, "sim_*/summary.txt"))

    for filepath in summary_files:
        data = {}
        sim_id = os.path.basename(os.path.dirname(filepath))
        data["batch"] = batch
        data["sim_id"] = sim_id
        with open(filepath) as f:
            for line in f:
                if "=" in line:
                    key, value = line.strip().split("=")
                    data[key.strip()] = value.strip()
        all_data.append(data)

df = pd.DataFrame(all_data)
df.to_csv(os.path.join(base_results_dir, "summary_table_all_batches.csv"), index=False)
