import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, "r") as h:
    noAstFaa = [x.replace("*", "") for x in h]

with open(output_file, "w") as h:
    for l in noAstFaa:
        h.write(f"{l}")