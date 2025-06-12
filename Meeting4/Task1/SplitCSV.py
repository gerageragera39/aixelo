import pandas as pd


def split_csv(path, name):

    df = pd.read_csv(path)

    df_co2 = df[df['mof_name'].str.endswith('_CO2')]
    df_h2o = df[df['mof_name'].str.endswith('_H2O')]

    df_co2.to_csv(f'{name}_CO2.csv', index=False)
    df_h2o.to_csv(f'{name}_H2O.csv', index=False)


if __name__ == '__main__':
    split_csv(r'C:\Users\samys\OneDrive\Рабочий стол\AXL_DATA\AXL\Meeting3\MyFingerprint\45print\45_fingerprints_grouped_avg_energy.csv', '45_fingerprints')
    split_csv(r'C:\Users\samys\OneDrive\Рабочий стол\AXL_DATA\AXL\Meeting3\MyFingerprint\120print\120_fingerprints_grouped_avg_energy.csv', '120_fingerprints')
