import numpy as np
import pandas as pd

import matplotlib as mpl
from matplotlib import colors 
import matplotlib.pyplot as plt

import seaborn as sns
import seaborn.objects as so

import sys, re


def extract_coordinates(row):
    positions = row['motif_coords']
    coordinates = re.findall(r'(\d+),(\d+)', positions)
    extracted_coords = []
    for start, end in coordinates:
        start = int(start)/row['seq_length']
        end = int(end)/row['seq_length']
        extracted_coords.append((start, row['seq_length']))
        extracted_coords.append((end, row['seq_length']))
        i = start + 0.005
        while i-0.005 < end:
            extracted_coords.append((i, row['seq_length']))
            i += 0.005
    return extracted_coords

def extract_signal(row):
    length = row['seq_length']
    extracted_coords = []
    start = int(row['signal_start']) / length
    end = int(row['signal_end']) / length
    extracted_coords.append((start, length))
    extracted_coords.append((end, length))
    i = start + 0.005
    while i-0.005 < end:
        extracted_coords.append((i, length))
        i += 0.005
    return extracted_coords


df = pd.read_table(sys.argv[1], header = None).dropna(axis=0)
config = sys.argv[1].split('.')[0]

coords = pd.DataFrame({'name': [], 'length': [], 'position': []})
for name, row in df.iterrows():
    if config == "_":
        coord = {'name': row[0], 'length': row[1], 'position': 0}
        coords = coords._append(coord, ignore_index=True)
    else:
        for i in range(2,len(row)):
            coord = {'name': row[0], 'length': row[1], 'position': row[i]}
            coords = coords._append(coord, ignore_index=True)
max_length = 1500



e_phenotype = coords[coords['name'].str.contains('effector|reduced_virulence|loss_of_pathogenicity')]
e_y_cys = e_phenotype.iloc[:, 1].to_numpy()
e_x_cys = e_phenotype.iloc[:, 2].to_numpy()

# rv_phenotype = coords[coords['name'].str.contains('reduced_virulence')]
# rv_y_cys = rv_phenotype.iloc[:, 1].to_numpy()
# rv_x_cys = rv_phenotype.iloc[:, 2].to_numpy()

# lop_phenotype = coords[coords['name'].str.contains('loss_of_pathogenicity')]
# lop_y_cys = lop_phenotype.iloc[:, 1].to_numpy()
# lop_x_cys = lop_phenotype.iloc[:, 2].to_numpy()

control_phenotype = coords[~coords['name'].str.contains('effector|reduced_virulence|loss_of_pathogenicity')]
control_y_cys = control_phenotype.iloc[:, 1].to_numpy()
control_x_cys = control_phenotype.iloc[:, 2].to_numpy()







e_list = e_phenotype['name'].drop_duplicates().tolist()
# rv_list = rv_phenotype['name'].drop_duplicates().tolist()
# lop_list = lop_phenotype['name'].drop_duplicates().tolist()
control_list = control_phenotype['name'].drop_duplicates().tolist()



cpp = pd.read_table(sys.argv[2], header = None).dropna(axis=0)
cpp_names = cpp.iloc[:, 0].tolist()

e_cpp = cpp[cpp[0].isin(e_list)]
e_y_cpp = e_cpp.iloc[:, 3].to_numpy()
e_x_cpp = e_cpp.iloc[:, 1].to_numpy()

# rv_cpp = cpp[cpp[0].isin(rv_list)]
# rv_y_cpp = rv_cpp.iloc[:, 3].to_numpy()
# rv_x_cpp = rv_cpp.iloc[:, 1].to_numpy()

# lop_cpp = cpp[cpp[0].isin(lop_list)]
# lop_y_cpp = lop_cpp.iloc[:, 3].to_numpy()
# lop_x_cpp = lop_cpp.iloc[:, 1].to_numpy()

control_cpp = cpp[cpp[0].isin(control_list)]
control_y_cpp = control_cpp.iloc[:, 3].to_numpy()
control_x_cpp = control_cpp.iloc[:, 1].to_numpy()



motifs = pd.read_table(sys.argv[3]).dropna(axis=0)

motifs['seq_length'] = motifs['seq'].apply(lambda x: len(x))
# motifs['cutsite'] = motifs[['signal_start', 'signal_end', 'kex2_start', 'kex2_end']].max(axis='columns')
# motifs['relative_cutsite'] = motifs['cutsite'] / motifs['seq_length']
motifs['relative_cutsite'] = motifs['kex2_end'] / motifs['seq_length']
motifs['relative_signal_peptide'] = motifs['signal_end'] / motifs['seq_length']

motifs['coordinates'] = motifs.apply(extract_coordinates, axis=1)
motifs['signal_coordinates'] = motifs.apply(extract_signal, axis=1)
sorted_motifs = motifs.sort_values(by=['seq_length'])

e_cutsites = motifs[motifs['name'].isin(e_list)]
e_y_cutsites = e_cutsites['seq_length'].to_numpy()
e_x_cutsites = e_cutsites['relative_cutsite'].to_numpy()
e_x_signal = [coord[0] for row in e_cutsites['signal_coordinates'] for coord in row]
e_y_signal = [coord[1] for row in e_cutsites['signal_coordinates'] for coord in row]
print(e_x_signal)
print(e_y_signal)
e_x_motifs = [coord[0] for row in e_cutsites['coordinates'] for coord in row]
e_y_motifs = [coord[1] for row in e_cutsites['coordinates'] for coord in row]


# rv_cutsites = motifs[motifs['name'].isin(rv_list)]
# rv_y_cutsites = rv_cutsites['seq_length'].to_numpy()
# rv_x_cutsites = rv_cutsites['relative_cutsite'].to_numpy()
# rv_x_signal = [coord[0] for row in rv_cutsites['signal_coordinates'] for coord in row]
# rv_y_signal = [coord[1] for row in rv_cutsites['signal_coordinates'] for coord in row]
# rv_x_motifs = [coord[0] for row in rv_cutsites['coordinates'] for coord in row]
# rv_y_motifs = [coord[1] for row in rv_cutsites['coordinates'] for coord in row]

# lop_cutsites = motifs[motifs['name'].isin(lop_list)]
# sorted_cutsites = lop_cutsites.sort_values(by=['seq_length']).drop(columns=['seq'])
# sorted_cutsites.to_csv(sys.argv[4] + ".intermediate.tsv", sep="\t")
# lop_y_cutsites = lop_cutsites['seq_length'].to_numpy()
# lop_x_cutsites = lop_cutsites['relative_cutsite'].to_numpy()
# lop_x_signal = [coord[0] for row in lop_cutsites['signal_coordinates'] for coord in row]
# lop_y_signal = [coord[1] for row in lop_cutsites['signal_coordinates'] for coord in row]
# lop_x_motifs = [coord[0] for row in lop_cutsites['coordinates'] for coord in row]
# lop_y_motifs = [coord[1] for row in lop_cutsites['coordinates'] for coord in row]

control_cutsites = motifs[motifs['name'].isin(control_list)]
control_y_cutsites = control_cutsites['seq_length'].to_numpy()
control_x_cutsites = control_cutsites['relative_cutsite'].to_numpy()
control_x_signal = [coord[0] for row in control_cutsites['signal_coordinates'] for coord in row]
control_y_signal = [coord[1] for row in control_cutsites['signal_coordinates'] for coord in row]
control_x_motifs = [coord[0] for row in control_cutsites['coordinates'] for coord in row]
control_y_motifs = [coord[1] for row in control_cutsites['coordinates'] for coord in row]

combined = pd.concat([e_cutsites, control_cutsites])
combined = combined.drop(columns=['relative_cutsite', 'relative_signal_peptide', 'signal_coordinates', 'seq', 'coordinates'])
combined['has_cpp'] = combined['name'].apply(lambda x: 1 if x in cpp_names else 0)
combined['has_signal'] = combined['signal_end'].apply(lambda x: 0 if x == 0 else 1)
combined['has_kex2'] = combined['kex2_end'].apply(lambda x: 0 if x == 0 else 1)
combined['num_cys'] = df.shape[1] - 2
combined.to_csv(sys.argv[4] + ".intermediate.tsv", sep="\t", index=False)


fig, axs = plt.subplots(2, figsize=(55, 25))

e_cys = axs[0].hist2d(e_x_cys, e_y_cys, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Blues', alpha=0.5)
e_cpp = axs[0].hist2d(e_x_cpp, e_y_cpp, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Reds', alpha=0.5)
e_cut = axs[0].hist2d(e_x_cutsites, e_y_cutsites, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Greys', alpha=0)
e_motifs = axs[0].hist2d(e_x_motifs, e_y_motifs, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Greens', alpha=0)
e_signal = axs[0].hist2d(e_x_signal, e_y_signal, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Purples', alpha=0.5)
axs[0].set_ylim([0, max_length])
axs[0].set_xlim([0, 1])
axs[0].set_title('effector_phenotype')
cbar1 = plt.colorbar(e_cys[3], ax=axs[0])
cbar1.set_label('Cys density')
cbar2 = plt.colorbar(e_cpp[3], ax=axs[0])
cbar2.set_label('CPP density')
cbar3 = plt.colorbar(e_cut[3], ax=axs[0])
cbar3.set_label('Cutsite density')
cbar4 = plt.colorbar(e_motifs[3], ax=axs[0])
cbar4.set_label('Motif density')
cbar5 = plt.colorbar(e_signal[3], ax=axs[0])
cbar5.set_label('Signal density')
axs[0].set_ylabel('Length')
axs[0].set_xlabel('Position relative to length')

# rv_cys = axs[0, 1].hist2d(rv_x_cys, rv_y_cys, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Blues', alpha=0.5)
# rv_cpp = axs[0, 1].hist2d(rv_x_cpp, rv_y_cpp, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Reds', alpha=0.5)
# rv_cut = axs[0, 1].hist2d(rv_x_cutsites, rv_y_cutsites, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Greys', alpha=0.5)
# rv_motifs = axs[0, 1].hist2d(rv_x_motifs, rv_y_motifs, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Greens', alpha=0.5)
# rv_signal = axs[0, 1].hist2d(rv_x_signal, rv_y_signal, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Purples', alpha=0.5)
# axs[0, 1].set_ylim([0, max_length])
# axs[0, 1].set_xlim([0, 1])
# axs[0, 1].set_title('reduced_virulence_phenotype')
# cbar6 = plt.colorbar(rv_cys[3], ax=axs[0, 1])
# cbar6.set_label('Cys density')
# cbar7 = plt.colorbar(rv_cpp[3], ax=axs[0, 1])
# cbar7.set_label('CPP density')
# cbar8 = plt.colorbar(rv_cut[3], ax=axs[0, 1])
# cbar8.set_label('Cutsite density')
# cbar9 = plt.colorbar(rv_motifs[3], ax=axs[0, 1])
# cbar9.set_label('Motif density')
# cbar10 = plt.colorbar(rv_signal[3], ax=axs[0, 1])
# cbar10.set_label('Signal density')
# axs[0, 1].set_ylabel('Length')
# axs[0, 1].set_xlabel('Position relative to length')


# lop_cys = axs[1, 0].hist2d(lop_x_cys, lop_y_cys, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmap='Blues', cmin = 1, alpha=0.5)
# lop_cpp = axs[1, 0].hist2d(lop_x_cpp, lop_y_cpp, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmap='Reds', cmin = 1, alpha=0.5)
# lop_cut = axs[1, 0].hist2d(lop_x_cutsites, lop_y_cutsites, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)],cmin = 1, cmap='Greys', alpha=0.5)
# lop_motifs = axs[1, 0].hist2d(lop_x_motifs, lop_y_motifs, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Greens', alpha=0.5)
# lop_signal = axs[1, 0].hist2d(lop_x_signal, lop_y_signal, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Purples', alpha=0.5)
# axs[1, 0].set_ylim([0, max_length])
# axs[1, 0].set_xlim([0, 1])
# axs[1, 0].set_title('loss_of_pathogenicity_phenotype')
# cbar11 = plt.colorbar(lop_cys[3], ax=axs[1, 0])
# cbar11.set_label('Cys density')
# cbar12 = plt.colorbar(lop_cpp[3], ax=axs[1, 0])
# cbar12.set_label('CPP density')
# cbar13 = plt.colorbar(lop_cut[3], ax=axs[1, 0])
# cbar13.set_label('Cutsite density')
# cbar14 = plt.colorbar(lop_motifs[3], ax=axs[1, 0])
# cbar14.set_label('Motif density')
# cbar15 = plt.colorbar(lop_signal[3], ax=axs[1, 0])
# cbar15.set_label('Signal density')
# axs[1, 0].set_ylabel('Length')
# axs[1, 0].set_xlabel('Position relative to length')


control_cys = axs[1].hist2d(control_x_cys, control_y_cys, bins = [np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Blues', alpha=0.5)
control_cpp = axs[1].hist2d(control_x_cpp, control_y_cpp, bins = [np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Reds', alpha=0.5)
control_cut = axs[1].hist2d(control_x_cutsites, control_y_cutsites, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Greys', alpha=0)
control_motifs = axs[1].hist2d(control_x_motifs, control_y_motifs, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Greens', alpha=0)
control_signal = axs[1].hist2d(control_x_signal, control_y_signal, bins=[np.arange(0, 1, 0.005), np.arange(0, max_length, 25)], cmin = 1, cmap='Purples', alpha=0.5)
axs[1].set_ylim([0, max_length])
axs[1].set_xlim([0, 1])
axs[1].set_title('control_phenotype')
cbar16 = plt.colorbar(control_cys[3], ax=axs[1])
cbar16.set_label('Cys density')
cbar17 = plt.colorbar(control_cpp[3], ax=axs[1])
cbar17.set_label('CPP density')
cbar18 = plt.colorbar(control_cut[3], ax=axs[1])
cbar18.set_label('Cutsite density')
cbar19 = plt.colorbar(control_motifs[3], ax=axs[1])
cbar19.set_label('Motif density')
cbar20 = plt.colorbar(control_signal[3], ax=axs[1])
cbar20.set_label('Signal density')
axs[1].set_ylabel('Length')
axs[1].set_xlabel('Position relative to length')

plt.savefig(sys.argv[4] + ".png")