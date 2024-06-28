import pandas as pd
import numpy as np
import sys

def find_kmer_positions(protein_sequence, kmer_length=4, required_aa={'D', 'R', 'E', 'K'}, threshold=3):

    positions = []
    sequence_length = len(protein_sequence)
    

    for i in range(sequence_length - kmer_length + 1):
        kmer = protein_sequence[i:i + kmer_length]

        count = sum(1 for aa in kmer if aa in required_aa)

        if count >= threshold:
            start = i + 1  
            end = i + kmer_length  
            positions.append(f"{start},{end}")

    if not positions:
        return np.nan

    return ';'.join(positions)


bed = pd.read_table(sys.argv[1], header=None, delimiter='\t').dropna(axis=0)
bed.columns = ['name', 'start', 'end']

intermediate = pd.read_table(sys.argv[2], delimiter='\t')

cpp_concatenated = pd.DataFrame({'name': [], 'cpp_coords': []})
gene_groups = bed.groupby(by=['name'])
for name, group in gene_groups:
    cpps = []
    group = group.sort_values(by=['start'])
    for i in range(len(group['start'])):
        coords = str(group['start'].values[i]) + "," + str(group['end'].values[i])
        cpps.append(coords)
    all_coords = ';'.join(cpps)
    cpp_concatenated = cpp_concatenated._append({'name': name, 'cpp_coords': all_coords}, ignore_index=True)

intermediate = intermediate.merge(cpp_concatenated, on=['name'], how='left').fillna('.')

DREK_concatenated = pd.read_table(sys.argv[3], delimiter='\t')
DREK_concatenated = DREK_concatenated.drop(['threshold', 'query_residues', 'motifs', 'ratios', 'seq'], axis=1)
DREK_concatenated.columns = ['name', 'DREK_coords']


# for name, row in intermediate.iterrows():
#     DREK_concatenated = DREK_concatenated._append({'name': row['name'], 'DREK_coords': find_kmer_positions(row['seq'])}, ignore_index=True)
    
intermediate = intermediate.merge(DREK_concatenated, on=['name'], how='left').fillna('.')

intermediate.to_csv("CPP_sites.positions.txt", sep="\t", index=False)