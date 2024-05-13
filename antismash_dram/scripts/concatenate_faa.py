from Bio import SeqIO
import csv
import sys

def concatenate_fasta(input_tsv, output_file):
    with open(output_file, 'w') as output_handle:
        with open(input_tsv, 'r') as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter='\t')
            for row in tsv_reader:
                fasta_file_path = row[0]
                with open(fasta_file_path, 'r') as fasta_handle:
                    for record in SeqIO.parse(fasta_handle, 'fasta'):
                        sequence = str(record.seq).replace('\n', '')
                        output_handle.write(f'>{record.id}\n{sequence}\n')

if __name__ == "__main__":
    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    concatenate_fasta(input_file_path, output_file_path)
