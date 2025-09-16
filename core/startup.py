'''
Entry point for the NCD corrections program.
This script initializes the program, prompts the user for input files, and starts the correction process.
'''
import core.CorrectionTools as CorrectionTools
import core.Corrector as Corrector
import numpy as np
import csv
def main():
    print("Welcome to the NCD corrections program!")
    
    '''Getting the user matrix file'''
    
    user_matrix_file = input("Please enter the path to the distance matrix file: ")
    #Creating an instance of CorrectionTools and correctorInstance to be used throughout program
    CorrectionToolsInstance = None
    correctorInstance = None
    try:
        
        # Creating the matrix from the input file
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
    user_output_format = input("TXT (1), CSV (2), Both (3):\n")

    #Creating the corrected distance matrix that the user had input
    print("Please choose a correction model:")
    user_model_input = input("Jukes Cantor (1), K80 (2), K81 (3), T92 (4), TN93 (5), All Models (6):\n")
    
    print("Please give a prefix to your corrected distance matrix:")
    user_prefix = input("Prefix: ")
    
    
    #Takes instance of the inputted matrix and runs the selected correction
    #NOTE: Implement feature that allows to have multiple corrections within a session
    try:
        
        def calculate_euclidean_distance(matrix1, matrix2):
            
            data_1 = [row[1:] for row in matrix1]
            
            data_2 = [row[1:] for row in matrix2]
            
            data_1 = np.array(data_1, dtype=float)
            data_2 = np.array(data_2, dtype=float)
            print(str(data_1.shape), str(data_2.shape))
            if data_1.shape != data_2.shape:
                print(data_1.shape, data_2.shape)
                raise ValueError("Matrices must have the same dimensions for Euclidean distance calculation.")
            
            distance = np.sqrt(np.sum((data_1 - data_2) ** 2))/(data_1.shape[0] ** 2 - data_1.shape[0]) #Modified to account for number of elements excluding diagonal for average distance
            
            return distance
        
        def correct_anomolies(matrix, originalDistances):
            small_counter = 0
            high_counter = 0
            one_counter = 0  # Counter for values equal to 1
            
            # Does not touch the diagonal on purpose
            for i in range(matrix.shape[0]):
                for j in range(matrix.shape[1]):
                    if i != j:
                        # Check for values equal to 1 FIRST (before any modifications)
                        if matrix[i][j] == 1.0:
                            print(f"1.0 at {i}, {j}")
                            one_counter += 1
                            matrix[i][j] = 0.85
                        
                        # Check for values less than or equal to 0.1
                        elif matrix[i][j] <= 0.1:
                            print(f"Less than 0.1 at {i}, {j}")
                            small_counter += 1
                            matrix[i][j] = 0.01
                        
                        # Check for large differences from original
                        elif abs(originalDistances[i][j] - matrix[i][j]) >= 0.20:
                            print(f"Original difference: {originalDistances[i][j] - matrix[i][j]}")
                            print(f"Difference too large, calculating average at {i}, {j}")
                            high_counter += 1
                            # Calculate average without abs() - remove if not needed
                            new_value = (matrix[i][j] + originalDistances[i][j]) / 2
                            matrix[i][j] = new_value
                            
                            # Check if the averaged value becomes 1.0 and correct it
                            if matrix[i][j] >= 0.98:
                                print(f"Averaged value became 1.0 at {i}, {j}, correcting to 0.85")
                                matrix[i][j] = 0.85
                                one_counter += 1
                    # Final safety check - catch any remaining values >= 0.99
                    if matrix[i][j] >= 0.98:
                        print(f"FINAL CATCH: Value >= 0.99 at {i}, {j}: {matrix[i][j]} -> 0.85")
                        matrix[i][j] = 0.85
                        one_counter += 1
            print(f"Amount of small values corrected: {small_counter}")
            print(f"Amount of Very Large Differences: {high_counter}")
            print(f"Amount of 1.0 values corrected: {one_counter}")
            return matrix
        
        def write_corrected_matrix_to_file(matrix, filename, formatnum, inputMatrix):
            # FINAL ANOMALY CORRECTION - Extract numeric data and apply correction
            labels = [row[0] for row in matrix]
            numeric_data = np.array([row[1:] for row in matrix])
            original_numeric_data = np.array([row[1:] for row in inputMatrix])
            
            # Apply final correction to ensure no 1.0 values
            final_corrected_data = correct_anomolies(matrix=numeric_data, originalDistances=original_numeric_data)
            
            # Reconstruct matrix with labels
            final_matrix = [[label] + list(row) for label, row in zip(labels, final_corrected_data)]
            
            euclidean_distance_file = open(filename.replace(".txt", "_euclidean_distance.txt"), "w")
            if formatnum == "1":
                correctedMatrixFile = open(filename, "w")
            elif formatnum == "2":
                correctedMatrixFile = open(filename.replace(".txt", ".csv"), "w")
            elif formatnum == "3":
                correctedMatrixFile = open(filename, "w")
                correctedMatrixFile2 = open(filename.replace(".txt", ".csv"), "w")
            else:
                print("Invalid format number")
                return

            correctedMatrixFile.write(str(len(final_matrix)) + "\n")  # Write the number of taxa

            for row in final_matrix:
                label = f"{row[0]:<10} "  # Left-align label, width 10
                if formatnum == "1":
                    distances = " ".join(f"{float(x):.6f}" for x in row[1:])
                    correctedMatrixFile.write(f"{label}{distances}\n")
                    
                elif formatnum == "2":
                    distances = ",".join(f"{float(x):.6f}" for x in row[1:])
                    correctedMatrixFile.write(f"{label}{distances}\n")
                    
                elif formatnum == "3":
                    distances = " ".join(f"{float(x):.6f}" for x in row[1:])  
                    correctedMatrixFile.write(f"{label}{distances}\n")
                    
                    distances2 = ",".join(f"{float(x):.6f}" for x in row[1:])
                    correctedMatrixFile2.write(f"{label}{distances2}\n")
                    
                else:
                    print("Invalid format number")
            
            euclidean_distance_file.write(f"The Euclidean distance between the two matrices is: {calculate_euclidean_distance(inputMatrix, final_matrix)}")
            euclidean_distance_file.close()
            
            if formatnum == "1":
                print(f"Matrix written to {filename}")
            elif formatnum == "2":
                print(f"Matrix written to {filename.replace('.txt', '.csv')}")
            elif formatnum == "3":
                print(f"Matrix written to {filename}")
                print(f"Matrix written to {filename.replace('.txt', '.csv')}")
                correctedMatrixFile2.close()
            correctedMatrixFile.close()


        model_map = {
            "1": (correctorInstance.jukesCantor, f"correction_output/{user_prefix}correctedJukesMatrixOptimized.txt"),
            "2": (correctorInstance.K80,   f"correction_output/{user_prefix}correctedK80MatrixOptimized.txt"),
            "3": (correctorInstance.K81,   f"correction_output/{user_prefix}correctedK81MatrixOptimized.txt"),
            "4": (correctorInstance.T92,   f"correction_output/{user_prefix}correctedT92MatrixOptimized.txt"),
            "5": (correctorInstance.TN93,  f"correction_output/{user_prefix}correctedTN93MatrixOptimized.txt"),
        }

        if user_model_input in model_map:
            model_func, out_file = model_map[user_model_input]
            matrix = model_func()
            normalized = CorrectionToolsInstance.normalize(inputMatrix=matrix)
            write_corrected_matrix_to_file(normalized, out_file, user_output_format, inputMatrix)
        elif user_model_input == "6":
            for model_func, out_file in model_map.values():
                matrix = model_func()
                normalized = CorrectionToolsInstance.normalize(inputMatrix=matrix)
                write_corrected_matrix_to_file(normalized, out_file, user_output_format, inputMatrix)
        else:
            print("Invalid input. Please try a number 1 - 6")
    except Exception as e:
        print(f"An error occurred during correction: {e}")


if __name__ == "__main__":
    main()