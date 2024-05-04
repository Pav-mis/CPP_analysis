import sys, os, re
from Bio import SeqIO

table = open(sys.argv[2], "w+") 
with open(sys.argv[1], "r") as fasta:
    for record in SeqIO.parse(fasta, "fasta"):
        row = str(record.id)+'\t'+str(record.seq)+'\n'
        table.write(row)
table.close