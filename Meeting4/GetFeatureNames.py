import pandas as pd

sample_path = r"C:\Users\samys\OneDrive\Рабочий стол\AXL_DATA\OneDrive_1_24.05.2025\stoich45_fingerprints.csv"  # ЗАМЕНИ НА СВОЙ ПУТЬ
reference_df = pd.read_csv(sample_path)
expected_columns = reference_df.columns.tolist()

sample_path2 = r"C:\Users\samys\OneDrive\Рабочий стол\AXL_DATA\AXL\Meeting4\generated_QMOF_45_fingerprints.csv"  # ЗАМЕНИ НА СВОЙ ПУТЬ
reference_df2 = pd.read_csv(sample_path2)
expected_columns2 = reference_df2.columns.tolist()

print(', '.join(f'{name} - {index}' for name, index in zip(expected_columns, range(len(expected_columns)))))
print(', '.join(f'{name} - {index}' for name, index in zip(expected_columns2, range(len(expected_columns2)))))
