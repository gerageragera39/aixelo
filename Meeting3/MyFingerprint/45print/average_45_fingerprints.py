import pandas as pd
import re

df = pd.read_csv("45_fingerprints_with_energy.csv")

def extract_base_mof_name(filename):
    return re.sub(r'(_random_\d+|_\d+)?_relaxed\.cif$', '', filename)

df['mof_name'] = df['filename'].map(extract_base_mof_name)

fingerprint_cols = [col for col in df.columns if col not in ['energy', 'filename', 'mof_name']]

agg_dict = {col: 'first' for col in fingerprint_cols}
agg_dict['energy'] = 'mean'

df_grouped = df.groupby('mof_name').agg(agg_dict).reset_index()

df_grouped.to_csv("45_fingerprints_grouped_avg_energy.csv", index=False)
print("Сохранено в 45_fingerprints_grouped_avg_energy.csv")