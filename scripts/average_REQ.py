import re
import os
def extract_req_values(verbose_output_file):
    req_values = []
    with open(verbose_output_file, 'r') as f:
        for line in f:
            # Look for lines with Re= pattern - improved regex for decimal numbers
            
            match = re.search(r'Re=(\d+.\d+)', line)
            if match:
                print(f"Found: {match.group(1)}")
                req_values.append(float(match.group(1)))
    return req_values

def calculate_average_req(req_values):
    return sum(req_values) / len(req_values) if req_values else 0

# Usage
input_dir = '/Users/damianmoreno/Downloads/NCD_Corrections_Program/REQ/src'
# Get all entries in the directory
all_entries = os.listdir(input_dir)
files_in_directory = [entry for entry in all_entries if os.path.isfile(os.path.join(input_dir, entry))]

output_dir = "REQ_AVG/Arabidopsis"
os.makedirs(output_dir, exist_ok=True)

for file in files_in_directory:
    if file.startswith('arabidopsis') and file.endswith('req_verbose_output.txt'):
        full_file_path = os.path.join(input_dir, file)
        req_values = extract_req_values(full_file_path)
        average_req = calculate_average_req(req_values)
        with open(os.path.join(output_dir, f"{file.replace('_req_verbose_output.txt', '_average_req_output.txt')}"), "w") as output_file:
            output_file.write(f"Average REQ: {average_req}\n")
            output_file.write(f"Number of values found: {len(req_values)}\n")
            
        output_file.close()
print("Averages Calculated Successfully")