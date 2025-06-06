import os
import numpy as np
import pandas as pd
from pymatgen.core import Structure, Element

ELEMENTS = [Element.from_Z(z).symbol for z in range(1, 104)]
element_names = [f'frac_{el}' for el in ELEMENTS]
chem_names = [
    "mean_atomic_weight","mean_column","mean_row","range_atomic_number","mean_atomic_number",
    "range_atomic_radius","mean_atomic_radius","range_electronegativity","mean_electronegativity",
    "avg_s_valence","avg_p_valence","avg_d_valence","avg_f_valence",
    "frac_s_valence","frac_p_valence","frac_d_valence","frac_f_valence"
]
final_columns = element_names + chem_names

def get_elemental_fractions(element_list):
    total = len(element_list)
    counts = {el: 0 for el in ELEMENTS}
    for el in element_list:
        if el in counts:
            counts[el] += 1
    fractions = [counts[el] / total if total > 0 else 0 for el in ELEMENTS]
    return fractions

def get_additional_chem_features(elements):
    Zs = np.array([Element(el).Z for el in elements])
    weights = np.array([Element(el).atomic_mass for el in elements])
    cols = np.array([Element(el).group if Element(el).group is not None else 0 for el in elements])
    rows = np.array([Element(el).row if getattr(Element(el), "row", None) is not None else 0 for el in elements])
    radii = np.array([Element(el).atomic_radius if Element(el).atomic_radius is not None else 0 for el in elements])
    enegs = np.array([Element(el).X if Element(el).X is not None else 0 for el in elements])

    s_v, p_v, d_v, f_v = [], [], [], []
    for el in elements:
        e = Element(el)
        try:
            s_v.append(e.full_electronic_structure.count(('s', 1)))
            p_v.append(e.full_electronic_structure.count(('p', 1)))
            d_v.append(e.full_electronic_structure.count(('d', 1)))
            f_v.append(e.full_electronic_structure.count(('f', 1)))
        except Exception:
            s_v.append(0)
            p_v.append(0)
            d_v.append(0)
            f_v.append(0)
    s_v, p_v, d_v, f_v = np.array(s_v), np.array(p_v), np.array(d_v), np.array(f_v)

    denom = np.sum(s_v + p_v + d_v + f_v) + 1e-10
    features = [
        np.mean(weights) if len(weights) else 0,
        np.mean(cols) if len(cols) else 0,
        np.mean(rows) if len(rows) else 0,
        np.ptp(Zs) if len(Zs) else 0,      # Range
        np.mean(Zs) if len(Zs) else 0,
        np.ptp(radii) if len(radii) else 0,
        np.mean(radii) if len(radii) else 0,
        np.ptp(enegs) if len(enegs) else 0,
        np.mean(enegs) if len(enegs) else 0,
        np.mean(s_v) if len(s_v) else 0,
        np.mean(p_v) if len(p_v) else 0,
        np.mean(d_v) if len(d_v) else 0,
        np.mean(f_v) if len(f_v) else 0,
        np.sum(s_v)/denom,
        np.sum(p_v)/denom,
        np.sum(d_v)/denom,
        np.sum(f_v)/denom,
    ]
    return features

root_dir = r'C:\Users\samys\OneDrive\Рабочий стол\AXL_DATA\AXL\cifs\pristine_H2O\pristine_H2O'

rows = []
mof_names = []


for mof_folder in sorted(os.listdir(root_dir)):
    folder_path = os.path.join(root_dir, mof_folder)
    if not os.path.isdir(folder_path):
        continue
    cif_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.cif')]
    all_elements = []
    for cif in cif_files:
        try:
            struct = Structure.from_file(cif)
            all_elements.extend([site.specie.symbol for site in struct.sites])
        except Exception as e:
            print(f"Ошибка чтения файла {cif}: {e}")
    if not all_elements:
        print(f"Нет атомов в {mof_folder} — пропускаем")
        continue
    elemental = get_elemental_fractions(all_elements)
    chem = get_additional_chem_features(all_elements)
    fingerprint = elemental + chem
    rows.append(fingerprint)
    mof_names.append(mof_folder)

df = pd.DataFrame(rows, columns=final_columns, index=mof_names)
df.index.name = "MOF"
df.to_csv("MOF_fingerprints.csv")
print("MOF_fingerprints.csv")