import pandas as pd
import sys

df = pd.read_table(sys.argv[1], header=None, delimiter='\t').dropna(axis=0)
df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "qcovs"]
ident_filtered = df[df['pident'] > 90.00]
cov_filtered = ident_filtered[ident_filtered['qcovs'] > 90.00]
bed = ident_filtered[['sseqid', 'sstart', 'send']]
print(bed)

bed.to_csv("CPP_sites.bed", sep="\t", header=False, index=False)