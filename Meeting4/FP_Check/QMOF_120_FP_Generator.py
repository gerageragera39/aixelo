import os
import pandas as pd
import numpy as np
from pymatgen.core import Structure, Element
from tqdm import tqdm

# Загрузка соответствия ID → имя MOF
name_df = pd.read_csv(r"C:\users\samys\Downloads\13147324\qmof_database\qmof_database\qmof.csv")
id_to_name = dict(zip(name_df["qmof_id"], name_df["name"]))

# Путь к папке с CIF-файлами
folder_path = r"C:\users\samys\Downloads\13147324\qmof_database\qmof_database\relaxed_structures\relaxed_structures"
cif_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".cif")]

# Список всех элементов
ELEMENTS = [Element.from_Z(z).symbol for z in range(1, 104)]
element_names = [f'frac_{el}' for el in ELEMENTS]
chem_names = [
    "mean_atomic_weight", "mean_column", "mean_row", "range_atomic_number", "mean_atomic_number",
    "range_atomic_radius", "mean_atomic_radius", "range_electronegativity", "mean_electronegativity",
    "avg_s_valence", "avg_p_valence", "avg_d_valence", "avg_f_valence",
    "frac_s_valence", "frac_p_valence", "frac_d_valence", "frac_f_valence"
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
        np.mean(weights),
        np.mean(cols),
        np.mean(rows),
        np.ptp(Zs),
        np.mean(Zs),
        np.ptp(radii),
        np.mean(radii),
        np.ptp(enegs),
        np.mean(enegs),
        np.mean(s_v),
        np.mean(p_v),
        np.mean(d_v),
        np.mean(f_v),
        np.sum(s_v)/denom,
        np.sum(p_v)/denom,
        np.sum(d_v)/denom,
        np.sum(f_v)/denom,
    ]
    return features

fingerprints = []
mof_names = []

# Обработка каждого CIF-файла
for cif_path in tqdm(cif_files):
    try:
        struct = Structure.from_file(cif_path)
        elements = [site.specie.symbol for site in struct.sites]
        elemental = get_elemental_fractions(elements)
        chem = get_additional_chem_features(elements)
        fingerprint = elemental + chem
    except Exception as e:
        print(f"Ошибка при обработке {cif_path}: {e}")
        fingerprint = [np.nan] * len(final_columns)

    qmof_id = os.path.splitext(os.path.basename(cif_path))[0]
    mof_name = id_to_name.get(qmof_id, qmof_id)
    mof_names.append(mof_name)
    fingerprints.append(fingerprint)

# Сбор в DataFrame
df = pd.DataFrame(fingerprints, columns=final_columns)
df.insert(0, "MOF", mof_names)
df.to_csv("generated_QMOF_120_fingerprints.csv", index=False)
