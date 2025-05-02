import os
import glob
import pandas as pd #type: ignore

batch_name = "test"
base_results_dir = os.path.expanduser("~/results")
batch_path = os.path.join(base_results_dir, batch_name)

summary_files = glob.glob(os.path.join(batch_path, "sim_*/summary.txt"))

all_summaries = []
for filepath in summary_files:
    data = {}
    sim_id = os.path.basename(os.path.dirname(filepath))
    data["sim_id"] = sim_id
    with open(filepath) as f:
        for line in f:
            if "=" in line:
                key, value = line.strip().split("=")
                data[key.strip()] = value.strip()
    all_summaries.append(data)

df = pd.DataFrame(all_summaries)
output_csv = os.path.join(batch_path, "summary_table.csv")
df.to_csv(output_csv, index=False)
