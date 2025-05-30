import os
import sys
import pandas as pd # type: ignore

if len(sys.argv) < 2:
    sys.exit(1)

# When using this command, you have to add the folder name after the command : "python merge_summary.py Folder_name"
# If you don't, the script will start merging all the folder of the "results" directory
folder_name = sys.argv[1]
base_results_dir = os.path.expanduser(f"~/results/{folder_name}")

# Search for "batch_" subfolders
batch_dirs = sorted([d for d in os.listdir(base_results_dir) if d.startswith("batch_")])

all_df = []

# Get the summary table inside of each subfolder 
for batch in batch_dirs:
    path = os.path.join(base_results_dir, batch, "summary_table.csv")
    if os.path.exists(path):
        if os.path.getsize(path) > 0:
            df = pd.read_csv(path)
            df["batch"] = batch
            all_df.append(df)
        else:
            print(f"Empty file, please check the summary at this location : {path}")
    else:
        pass

# Merging all summary table in order into one summary and put it into csv format
if all_df:
    merged = pd.concat(all_df, ignore_index=True)
    output_file = os.path.join(base_results_dir, "summary_table_all_batches.csv")
    merged.to_csv(output_file, index=False)
    print(f"Merge complete : {len(merged)} new lines in {output_file}")
else:
    pass
