import pandas as pd
import re

df = pd.read_csv("120_fingerprints_with_energy.csv", index_col=0)

def extract_base_mof_name(filename):
    return re.sub(r'(_random_\d+|_\d+)?_relaxed\.cif$', '', filename)

df['mof_name'] = df.index.map(extract_base_mof_name)

fingerprint_cols = [col for col in df.columns if col not in ['energy', 'mof_name']]
agg_dict = {col: 'first' for col in fingerprint_cols}
agg_dict['energy'] = 'mean'

df_grouped = df.groupby('mof_name').agg(agg_dict)
df_grouped.index.name = 'mof_name'

df_grouped.to_csv("120_fingerprints_grouped_avg_energy.csv")
print("Сохранено в 120_fingerprints_grouped_avg_energy.csv")