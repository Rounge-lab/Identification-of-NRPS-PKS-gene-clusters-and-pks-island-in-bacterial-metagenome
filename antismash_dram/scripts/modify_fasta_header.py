from Bio import SeqIO
import sys
import os

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

path = os.path.abspath(input_fasta)
string_to_add = os.path.splitext(path.rsplit("/")[-3])[0] + "_" + os.path.splitext(path.rsplit("/")[-2])[0] + "_" 

fasta_sequences = list(SeqIO.parse(input_fasta, "fasta"))

for record in fasta_sequences:
    record.id = string_to_add + record.id

with open(output_fasta, "w") as output_handle:
    SeqIO.write(fasta_sequences, output_handle, "fasta")