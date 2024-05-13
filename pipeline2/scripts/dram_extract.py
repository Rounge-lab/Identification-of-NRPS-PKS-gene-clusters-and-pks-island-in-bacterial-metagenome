from Bio import SeqIO
import pandas as pd
import os
import sys
import csv
import re

def read_dram_tsv(tsv_file):
    records = []
    with open(tsv_file, "r") as file:
        tsv_file = csv.reader(file, delimiter="\t")
        header = file.readline().strip().split('\t')
        for line in tsv_file:
            name = line[0].split("_")
            col = line[0]
            # gene_id
            gene_id_1 = re.compile(r"ctg\d+_\d+")
            gene_id_2 = re.compile(r"allorf_\d+_\d+")
            gene_id_match_1 = gene_id_1.search(col)
            gene_id_match_2 = gene_id_2.search(col)
            if gene_id_match_1:
                gene_id = gene_id_match_1.group(0)
            elif gene_id_match_2:
                gene_id = gene_id_match_2.group(0)
            # genome
            genome_1 = re.compile(r"_(S-\d*_\w*_\d*_sub)_")
            genome_2 = re.compile(r"_(S-\d*_\w*_\d*)_ctg")
            genome_3 = re.compile(r"_(S-\d*_\w*_\d*)_allorf")
            genome_match_1 = genome_1.search(col)
            genome_match_2 = genome_2.search(col)
            genome_match_3 = genome_3.search(col)
            if genome_match_1:
                genome = genome_match_1.group(1)
            elif genome_match_2:
                genome = genome_match_2.group(1)
            elif genome_match_3:
                genome = genome_match_3.group(1)
            # uniref_annotation
            match_anno = re.match(r"^UniRef90_\w*\s*(.*?)\s*n=\d+", line[header.index("uniref_hit")+1])
            if line[header.index("uniref_hit")+1] != "":
                annotation = match_anno.group(1)
            else:
                annotation = ""
            record = {
                "Genome": genome,
                "MAG": "_".join(name[3:4]),
                "Gene_ID": gene_id,
                "Uniref_ID": line[header.index("uniref_id")+1],
                "Uniref_Annotation": annotation,
                "Uniref_Taxonomy": line[header.index("uniref_taxonomy")+1],
                "Uniref_Identity": line[header.index("uniref_identity")+1],
                "Uniref_eVal": line[header.index("uniref_eVal")+1]
                # "Kegg_Annotation": line[4],
                # "Peptidase_Annotation": line[14],
                # "Pfam_Annotation": line[19],
                # "Cazy_Annotation": line[21]
            }
            records.append(record)
    return records

def main():
    input_file = sys.argv[1] # dram_tsv_output
    output_file = sys.argv[2] # dram_annotation_result

    dram = []
    dram.extend(read_dram_tsv(input_file))
    dram = pd.DataFrame(dram)

    dram.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()
