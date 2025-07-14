def convert_txt_to_csv_file(input_filepath, output_filepath):
    """
    Converts a text file with space-separated values into a CSV file.

    Args:
        input_filepath (str): The path to the input .txt file.
        output_filepath (str): The path where the output .csv file will be saved.
    """
    try:
        with open(input_filepath, 'r') as infile:
            file_content = infile.read()

        csv_lines = []
        for line in file_content.strip().split('\n'):
            # Skip empty lines
            if not line.strip():
                continue
            # Split each line by any whitespace and join with a comma
            csv_line = ','.join(line.split())
            csv_lines.append(csv_line)

        with open(output_filepath, 'w') as outfile:
            outfile.write('\n'.join(csv_lines))

        print(f"Successfully converted '{input_filepath}' to '{output_filepath}'")

    except FileNotFoundError:
        print(f"Error: The file '{input_filepath}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
        
convert_txt_to_csv_file('correction_output/correctedTN93MatrixTest.txt', 'csv_correction_output/correctedTN93MatrixTest.csv')