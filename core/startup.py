'''
Entry point for the NCD corrections program.
This script initializes the program, prompts the user for input files, and starts the correction process.
'''
import core.CorrectionTools as CorrectionTools
import core.Corrector as Corrector
import csv
def main():
    print("Welcome to the NCD corrections program!")
    
    '''Getting the user matrix file'''
    
    user_matrix_file = input("Please enter the path to the distance matrix file: ")
    #Creating an instance of CorrectionTools and correctorInstance to be used throughout program
    CorrectionToolsInstance = None
    correctorInstance = None
    try:
        
        # Creating the matrix form the input file
        CorrectionToolsInstance = CorrectionTools.CorrectionTools(user_matrix_file)
        inputMatrix = CorrectionToolsInstance.createMatrix()
            
    except FileNotFoundError:
        print("Error: File not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
        
    
    '''Getting user FASTA file'''
    
    user_FASTA_file = input("Please enter the path to the FASTA file: ")
    try:
        correctorInstance = Corrector.Corrector(inputMatrix, user_FASTA_file)
        #CorrectionToolsInstance 
    except FileNotFoundError:
        print("Error: File not found. Please enter a valid path.")
    except Exception as e:
        print(f"An error occurred: {e}")

    print("Please choose an output format for the corrected distance matrix:")
    user_output_format = input("TXT (1), CSV (2)\n")

    #Creating the corrected distance matrix that the user had input
    print("Please choose a correction model:")
    user_model_input = input("Jukes Cantor (1), K80 (2), K81 (3), T92 (4), TN93 (5), All Models (6):\n")
    
    
    #Takes instance of the inputted matrix and runs the selected correction
    #NOTE: Implement feature that allows to have multiple corrections within a session
    try:
        
        def write_corrected_matrix_to_file(matrix, filename, formatnum):
            if formatnum == "1":
                with open(filename, "w") as correctedMatrixFile:
                    correctedMatrixFile.write(str(len(matrix)) + "\n")  # Write the number of rows at top
                    for row in matrix:
                        correctedMatrixFile.write(" ".join(str(x) for x in row))
                        correctedMatrixFile.write("\n")
                print(f"Matrix written to {filename}")
            elif formatnum == "2":
                with open(filename.replace(".txt", ".csv"), "w") as correctedMatrixFile:
                    writer = csv.writer(correctedMatrixFile)
                    writer.writerows(matrix)
                print(f"Matrix written to {filename.replace(".txt", ".csv")}")

        model_map = {
            "1": (correctorInstance.jukesCantor, "correction_output/correctedJukesMatrixOptimized.txt"),
            "2": (correctorInstance.K80,   "correction_output/correctedK80MatrixOptimized.txt"),
            "3": (correctorInstance.K81,   "correction_output/correctedK81MatrixOptimized.txt"),
            "4": (correctorInstance.T92,   "correction_output/correctedT92MatrixOptimized.txt"),
            "5": (correctorInstance.TN93,  "correction_output/correctedTN93MatrixOptimized.txt"),
        }

        if user_model_input in model_map:
            model_func, out_file = model_map[user_model_input]
            matrix = model_func()
            normalized = CorrectionToolsInstance.normalize(inputMatrix=matrix)
            write_corrected_matrix_to_file(normalized, out_file, user_output_format)
        elif user_model_input == "6":
            for model_func, out_file in model_map.values():
                matrix = model_func()
                normalized = CorrectionToolsInstance.normalize(inputMatrix=matrix)
                write_corrected_matrix_to_file(normalized, out_file, user_output_format)
        else:
            print("Invalid input. Please try a number 1 - 6")
    except Exception as e:
        print(f"An error occurred during correction: {e}")


if __name__ == "__main__":
    main()