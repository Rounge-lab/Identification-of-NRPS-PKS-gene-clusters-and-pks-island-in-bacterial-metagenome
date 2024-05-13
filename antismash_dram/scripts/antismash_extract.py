from Bio import SeqIO
import pandas as pd
import os
import sys
import csv
import re

def read_anti_cds(gbk_file):
    records = []
    path = os.path.abspath(gbk_file)
    for seq_record in SeqIO.parse(gbk_file, "genbank"):
        for seq_feature in seq_record.features:
            if seq_feature.type == "CDS":
                assert len(seq_feature.qualifiers["translation"]) == 1
                record = {
                    "MAG": os.path.splitext(path.rsplit("/")[-3])[0],
                    "Genome": os.path.splitext(path.rsplit("/")[-2])[0],
                    "Region": os.path.splitext(path.rsplit("/")[-1])[0],
                    "Gene_ID": seq_feature.qualifiers["locus_tag"][0],
                    "Gene_Location_Start": seq_feature.location.start,
                    "Gene_Location_End": seq_feature.location.end,
                    "Sequence": seq_feature.qualifiers["translation"][0],
                    "NRPS_PKS_Modules": None
                }
                records.append(record)
                if "NRPS_PKS" in seq_feature.qualifiers:
                    nrps_pks_list = seq_feature.qualifiers["NRPS_PKS"]
                    nrps_pks = []
                    for element in nrps_pks_list:
                        match = re.match(r"Domain: ([\w\d_-]+)", element)
                        if match:
                            nrps_pks.append(match.group(1))
                    record["NRPS_PKS_Modules"] = nrps_pks
    return records

def read_anti_protocluster(gbk_file):
    records = []
    path = os.path.abspath(gbk_file)
    with open(gbk_file, 'r') as file:
        lines = file.readlines()
        for line in lines[1:20]:
            match_start = re.match(r"^\s*Orig\. start\s*::\s*(\d+)", line)
            match_end = re.match(r"^\s*Orig\. end\s*::\s*(\d+)", line)
            if match_start:
                start = match_start.group(1)
            if match_end:
                end = match_end.group(1)
    for seq_record in SeqIO.parse(gbk_file, "genbank"):
        for seq_feature in seq_record.features:
            if seq_feature.type == "protocluster":
                record = {
                    "MAG": os.path.splitext(path.rsplit("/")[-3])[0],
                    "Genome": os.path.splitext(path.rsplit("/")[-2])[0],
                    "Region": os.path.splitext(path.rsplit("/")[-1])[0],
                    "Original_Start": start,
                    "Original_End": end,
                    "Protocluster_Number": seq_feature.qualifiers["protocluster_number"][0],
                    "Protocluster_Category": seq_feature.qualifiers["category"][0],
                    "Protocluster_Product": seq_feature.qualifiers["product"][0],
                    "Protocluster_Location_Start": seq_feature.location.start,
                    "Protocluster_Location_End": seq_feature.location.end
                }
                records.append(record)
    return records

def read_anti_aSModule(gbk_file):
    records = []
    for seq_record in SeqIO.parse(gbk_file, "genbank"):
        for seq_feature in seq_record.features:
            if seq_feature.type == "aSModule":
                if "complete" in seq_feature.qualifiers:
                    record = {
                        "Gene_ID": seq_feature.qualifiers["locus_tags"][0],
                        "Amino_Acid": seq_feature.qualifiers["monomer_pairings"][0]
                        #"aSModule_Location_Start": seq_feature.location.start,
                        #"aSModule_Location_End": seq_feature.location.end
                    }
                    records.append(record)
    return records

def main():
    input_file = sys.argv[1] # antimash_gbk_output
    output_file_1 = sys.argv[2] # protocluster_level_antismash_result
    output_file_2 = sys.argv[3] # gene_level_antismash_result

    cds = []
    cds.extend(read_anti_cds(input_file))
    cds = pd.DataFrame(cds)

    protocluster = []
    protocluster.extend(read_anti_protocluster(input_file))
    protocluster = pd.DataFrame(protocluster)

    cluster_level = protocluster.copy()
    cluster_level["Gene_ID"] = [[] for _ in range(len(cluster_level))]

    for protocluster_index, protocluster_row in protocluster.iterrows():
        for cds_index, cds_row in cds.iterrows():
            if (cds_row["Gene_Location_Start"] >= protocluster_row["Protocluster_Location_Start"] and cds_row["Gene_Location_End"] <= protocluster_row["Protocluster_Location_End"]):
                cluster_level.at[protocluster_index, "Gene_ID"].append(cds_row["Gene_ID"])
    cluster_level.to_csv(output_file_1, sep="\t", index=False)

    asmodule = []
    asmodule.extend(read_anti_aSModule(input_file))
    asmodule = pd.DataFrame(asmodule)

    if not asmodule.empty:
        asmodule = asmodule.groupby("Gene_ID")["Amino_Acid"].agg(list).reset_index()
        gene_level = pd.merge(cds, asmodule, on="Gene_ID", how="outer")
        gene_level = pd.DataFrame(gene_level)
    else:
        cds["Amino_Acid"] = None
        gene_level = cds
        gene_level = pd.DataFrame(gene_level)
    gene_level.to_csv(output_file_2, sep='\t', index=False)

if __name__ == "__main__":
    main()
