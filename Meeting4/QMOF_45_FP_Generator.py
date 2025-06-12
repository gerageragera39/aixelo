import os
import pandas as pd
import numpy as np
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Element
from scipy.stats import gmean
from tqdm import tqdm

# Загрузка таблицы с именами MOF
name_df = pd.read_csv(r"C:\users\samys\Downloads\13147324\qmof_database\qmof_database\qmof.csv")
id_to_name = dict(zip(name_df["qmof_id"], name_df["name"]))

# Свойства элементов
property_mapping = {
    "atomic_num": "Z",
    "group": "group",
    "period": "row",
    "electronegativity": "X",
    "electron_affinity": "electron_affinity",
    "melting_point": "melting_point",
    "boiling_point": "boiling_point",
    "density": "density_of_solid",
    "ionization_energy": "ionization_energies"
}
stat_names = ["mean", "geometric_mean", "standard_deviation", "max", "min"]


def compute_stats(values):
    values = np.array(values, dtype=np.float64)
    values = values[~np.isnan(values)]
    if len(values) == 0:
        return [np.nan] * 5
    return [
        np.mean(values),
        gmean(values) if np.all(values > 0) else np.nan,
        np.std(values),
        np.max(values),
        np.min(values)
    ]


# 🔧 Путь к папке с CIF-файлами
folder_path = r"C:\users\samys\Downloads\13147324\qmof_database\qmof_database\relaxed_structures\relaxed_structures"
cif_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".cif")]

fingerprints = []
mof_names = []

# Обработка каждого CIF
for path in tqdm(cif_files):
    try:
        structure = Structure.from_file(path)
        elements = [site.specie for site in structure if isinstance(site.specie, Element)]
        stats = []
        for prop_key, attr in property_mapping.items():
            values = []
            for el in elements:
                try:
                    if attr == "ionization_energies":
                        val = el.ionization_energies[0] if el.ionization_energies else np.nan
                    else:
                        val = getattr(el, attr, np.nan)
                    values.append(val)
                except Exception:
                    values.append(np.nan)
            stats.extend(compute_stats(values))
    except Exception as e:
        print(f"Ошибка с файлом {path}: {e}")
        stats = [np.nan] * 45

    # Получаем имя файла без расширения → это qmof_id
    qmof_id = os.path.splitext(os.path.basename(path))[0]
    mof_name = id_to_name.get(qmof_id, qmof_id)  # Если не найдено — оставить как есть
    mof_names.append(mof_name)
    fingerprints.append(stats)

# Сбор результатов
columns = [f"{prop}_{stat}" for prop in property_mapping for stat in stat_names]
fingerprints_df = pd.DataFrame(fingerprints, columns=columns)
fingerprints_df.insert(0, "MOF", mof_names)
fingerprints_df.to_csv("generated_QMOF_45_fingerprints.csv", index=False)
