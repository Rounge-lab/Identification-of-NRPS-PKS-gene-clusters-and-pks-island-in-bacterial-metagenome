import sys
import pandas as pd
import numpy as np

def bedcov(input_file_path, output_file_path):
    with open(input_file_path, 'r') as file_list:
        with open(output_file_path, 'w') as output_file:
            for file_path in file_list:
                file_path = file_path.strip()
                try:
                    with open(file_path, 'r') as input_file:
                        df = pd.read_table(input_file, delim_whitespace = True, header = None)
                        bedcov = df[4]
                        df = pd.DataFrame(columns=bedcov)
                        df[file_path.split("/")[-1].split(".")[0]] = np.nan
                        df.to_csv(output_file, sep='\t', index=False)
                except FileNotFoundError:
                    print(f"File not found: {file_path}")
                except Exception as e:
                    print(f"Error reading file {file_path}: {str(e)}")

if __name__ == "__main__":
    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    bedcov(input_file_path, output_file_path)
