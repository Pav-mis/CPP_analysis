import pandas as pd
import numpy as np
import sys

intermediate = pd.read_csv(sys.argv[1], delimiter='\t')


for name, row in intermediate.iterrows():
    df = pd.DataFrame({
    'aa_number': list(range(1, len(row['seq'])+1)),
    'aa': list(row['seq']),
    'group': pd.Series([np.nan] * (len(row['seq'])))
    })
    df["group"] = df["group"].astype('string')

    print(row['name'])

    hydro = pd.read_csv(sys.argv[2] + '/' + row['name'] + '1.dat', delimiter='\t', header=None)
    hydro.columns = ['aa_number', 'hydropathy']
    print(hydro)

    df["group"] = df["group"].astype('string')
    df = pd.merge(df, hydro, on=['aa_number'], how='left')

    #add cys
    for f in range(1, len(row['seq'])+1):
        if df.at[f-1, "aa"] == 'C':
            df.at[f-1, "group"] = 'cysteine'

    #add signal peptide
    if row['signal_start'] != 0:
        for f in range(row['signal_start'], row['signal_end']):
            df.at[f-1, "group"] = 'signal_peptide'

    #add kex2
    if row['kex2_start'] != 0:
        for f in range(row['kex2_start'], row['kex2_end']):
            df.at[f-1, "group"] = 'KEX2'

    if row['cpp_coords'] != '.':
        cpp_sites = row['cpp_coords'].split(';')
        for site in cpp_sites:
            coords = site.split(',')
            for f in range(int(coords[0]), int(coords[1])):
                df.at[f-1, "group"] = 'CPP'
    
    if row['DREK_coords'] != '.':
        DREK_mers = row['DREK_coords'].split(';')
        for mer in DREK_mers:
            coords = mer.split(',')
            for f in range(int(coords[0]), int(coords[1])):
                if pd.isna(df.at[f-1, "group"]):
                    df.at[f-1, "group"] = 'DREK'
                elif df.at[f-1, "group"] != "CPP":
                    df.at[f-1, "group"] = 'DREK'
                else:
                    df.at[f-1, "group"] = 'DREK+CPP'
            

    df.to_csv(sys.argv[3]+"/" + row['name'] + ".csv", index=False)