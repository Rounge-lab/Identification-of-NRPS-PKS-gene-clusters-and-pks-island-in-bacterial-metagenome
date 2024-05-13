from Bio import SeqIO
import sys
import os

gbk_file = sys.argv[1]
faa_file = sys.argv[2]

if not os.path.exists(gbk_file):
    sys.exit(0)

input_handle  = open(gbk_file, "r")
output_handle = open(faa_file, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            output_handle.write(">%s from %s\n%s\n" % (
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_record.name,
                   seq_feature.qualifiers['translation'][0]))

output_handle.close()
input_handle.close()
