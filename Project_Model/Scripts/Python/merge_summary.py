import os
import pandas as pd # type: ignore

base_results_dir = os.path.expanduser("~/results/table_test")
batch_dirs = sorted([d for d in os.listdir(base_results_dir) if d.startswith("batch_table_ref_")])

all_df = []

for batch in batch_dirs:
    path = os.path.join(base_results_dir, batch, "summary_table.csv")
    if os.path.exists(path):
        if os.path.getsize(path) > 0:
            df = pd.read_csv(path)
            df["batch"] = batch
            all_df.append(df)
        else:
            print(f"⚠️ Fichier vide ignoré : {path}")
    else:
        pass

if all_df:
    merged = pd.concat(all_df, ignore_index=True)
    output_file = os.path.join(base_results_dir, "summary_table_all_batches.csv")
    merged.to_csv(output_file, index=False)
    print(f"Merge complete : {len(merged)} new lines in {output_file}")
else:
    pass
