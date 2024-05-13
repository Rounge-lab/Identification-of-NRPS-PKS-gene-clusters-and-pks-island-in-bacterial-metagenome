import sys
import os

input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

path = os.path.abspath(input_file_path)
mag = os.path.splitext(path.rsplit("/")[-4])[0]
genome = os.path.splitext(path.rsplit("/")[-3])[0]
contig = os.path.splitext(path.rsplit("/")[-1])[0]

with open(input_file_path, 'r') as input_file:
    # Open the output file for writing
    with open(output_file_path, 'w') as output_file:
        # Flag to indicate if we are currently processing relevant lines
        processing_hits = False
        hit_info = {}
        # Iterate through each line in the input file
        for line in input_file:
            # Check if the line starts with "Significant hits"
            if line.startswith("Significant hits"):
                # Set the flag to True, indicating we should start processing hits
                processing_hits = True
            elif processing_hits:
                # Check if the line is empty, indicating the end of significant hits
                if line.strip() == "":
                    # Stop processing hits once empty line is found
                    processing_hits = False
                else:
                    # Extract significant hits and map them to details
                    hit_number, bgc_info = line.strip().split('. ', 1)
                    bgc_accession, bgc_name = bgc_info.split('\t', 1)
                    hit_info[int(hit_number)] = {'bgc_accession': bgc_accession, 'bgc_name': bgc_name}
            # Check if the line starts with "Details"
            elif line.startswith("Details"):
                # Iterate through the lines in the details section for the current hit
                for line in input_file:
                    if line.startswith(">>"):
                        info_line = next(input_file).strip()
                        hit_number, bgc_accession = info_line.strip().split('. ', 1)
                    elif line.startswith("Type:"):
                        bgc_type = line.split(": ", 1)[1].strip()
                        hit_info[int(hit_number)]["bgc_type"] = bgc_type

        # Write the extracted information to the output file
        for hit_number, info in hit_info.items():
            output_file.write(f"{mag}\t{genome}\t{contig}\t{info['bgc_accession']}\t{info['bgc_name']}\t{info['bgc_type']}\n")
