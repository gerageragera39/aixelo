import csv

filename = r'C:\Users\samys\OneDrive\Рабочий стол\AXL_DATA\AXL\Meeting3\MyFingerprint\45print\45_fingerprints_grouped_avg_energy.csv'  # замените на имя вашего файла

with open(filename, newline='', encoding='utf-8') as csvfile:
    reader = csv.reader(csvfile)
    first_row = next(reader)  # читаем первую строку
    num_columns = len(first_row)

print(f"Количество столбцов в CSV: {num_columns}")
