import sys
import pandas as pd

def concatenate_files(input_file_path, output_file_path):
    with open(input_file_path, 'r') as file_list:
        with open(output_file_path, 'w') as output_file:
            for file_path in file_list:
                file_path = file_path.strip()
                try:
                    with open(file_path, 'r') as input_file:
                        output_file.write(input_file.read())
                except FileNotFoundError:
                    print(f"File not found: {file_path}")
                except Exception as e:
                    print(f"Error reading file {file_path}: {str(e)}")

if __name__ == "__main__":
    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    concatenate_files(input_file_path, output_file_path)
