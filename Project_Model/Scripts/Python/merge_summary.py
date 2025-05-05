import os
import pandas as pd # type: ignore

base_results_dir = os.path.expanduser("~/results")
batch_dirs = sorted([d for d in os.listdir(base_results_dir) if d.startswith("batch_")])

all_df = []

for batch in batch_dirs:
    path = os.path.join(base_results_dir, batch, "summary_table.csv")
    if os.path.exists(path):
        df = pd.read_csv(path)
        df["batch"] = batch  # Ajoute l’identifiant du batch
        all_df.append(df)
    else:
        print(f"⚠️ Aucun fichier trouvé pour {batch}")

if all_df:
    merged = pd.concat(all_df, ignore_index=True)
    output_file = os.path.join(base_results_dir, "summary_table_all_batches.csv")
    merged.to_csv(output_file, index=False)
    print(f"✅ Fusion réussie : {len(merged)} lignes écrites dans {output_file}")
else:
    print("❌ Aucun fichier summary_table.csv trouvé.")
