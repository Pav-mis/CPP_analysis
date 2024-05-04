import pandas as pd
import sys

ranked = pd.read_csv(sys.argv[1], delimiter='\t')
gff = pd.read_csv(sys.argv[2], delimiter='\t', header=None)
seqs = pd.read_csv(sys.argv[3], delimiter='\t', header=None)
seqs.columns = ['name', 'seq']


signal_P = gff[gff[1].str.contains("SignalP:5.0b")]
kex_2 = gff[gff[8].str.contains("kex2_cutsite")]
pfam = gff[gff[1].str.contains("Pfam-scan")]
signal_coords = signal_P[[0, 3, 4]]
signal_coords.columns = ['name', 'signal_start', 'signal_end']

seqs_and_signal = seqs.merge(signal_coords, on=['name'], how='left')
seqs_and_signal = seqs_and_signal.set_index('name').fillna(0)

kex2_after_signal = pd.DataFrame({'name': [], 'kex2_start': [], 'kex2_end': []})
print

kex_groups = kex_2.groupby(by=[0])
for name, group in kex_groups:
    name = seqs_and_signal.loc[name[0]]
    group = group.sort_values(by=[3])
    signal_end = name['signal_end']
    for index, row in group.iterrows():
        if row[3] > signal_end:
            kex2_after_signal = kex2_after_signal._append({'name': row[0], 'kex2_start': row[3], 'kex2_end': row[4]}, ignore_index=True)
            break
kex2_after_signal = kex2_after_signal.set_index('name')

seqs_and_signal_and_kex2 = seqs_and_signal.merge(kex2_after_signal, on=['name'], how='left').fillna(0)

pfam_concatenated = pd.DataFrame({'name': [], 'motif_coords': []})
pfam_groups = pfam.groupby(by=[0])
for name, group in pfam_groups:
    motifs = []
    group = group.sort_values(by=[3])
    for i in range(len(group[3])):
        coords = str(group[3].values[i]) + "," + str(group[4].values[i])
        motifs.append(coords)
    all_coords = ';'.join(motifs)
    pfam_concatenated = pfam_concatenated._append({'name': name[0], 'motif_coords': all_coords}, ignore_index=True)

seqs_and_signal_and_kex2_and_motifs = seqs_and_signal_and_kex2.merge(pfam_concatenated, on=['name'], how='left').fillna('.')
seqs_and_signal_and_kex2_and_motifs = seqs_and_signal_and_kex2_and_motifs.astype({'signal_start': int, 'signal_end': int, 'kex2_start': int, 'kex2_end': int})


seqs_and_signal_and_kex2_and_motifs.to_csv("intermediate.tsv", index=False, sep='\t')