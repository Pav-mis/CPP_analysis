import pandas as pd
import sys

bed = pd.read_table(sys.argv[1], header=None, delimiter='\t').dropna(axis=0)
bed.columns = ['name', 'start', 'end']
cys = pd.read_table(sys.argv[2], header=None, delimiter='\t').dropna(axis=0)
cys.columns = ['name', 'length', 'cys_count', 'config', 'distances']
cys = cys.drop(columns=['cys_count', 'config', 'distances'])
cys_dict = cys.set_index('name').to_dict('dict')
lengths = cys_dict['length']
bed['length'] = bed["name"].map(lengths)
bed['start'] = bed['start']/bed['length']
bed['end'] = bed['end']/bed['length']
print(bed)

bed.to_csv("CPP_sites.positions.txt", sep="\t", header=False, index=False)
