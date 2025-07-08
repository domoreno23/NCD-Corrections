'''
Entry point for the NCD corrections program.
This script initializes the program, prompts the user for input files, and starts the correction process.
'''
import core.CorrectionTools as CorrectionTools
import core.Corrector as Corrector

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

    #Creating the corrected distance matrix that the user had input
    print("Please choose a correction model:")
    user_model_input = input("Jukes Cantor (1), K80 (2), K81 (3), T92 (4), TN93 (5):\n")
    
    
    #Takes instance of the inputted matrix and runs the selected correction
    #NOTE: Implement feature that allows to have multiple corrections within a session
    try:
        if user_model_input == "1":
            jukesMatrix = correctorInstance.jukesCantor()
            normalizedJukesMatrix = CorrectionToolsInstance.normalize(inputMatrix= jukesMatrix)
            #Writing the corrected Matrix to an output file
            correctedMatrixFile = open("correction_output/correctedJukesMatrix.txt", "w")
            for i in range(len(normalizedJukesMatrix)):
                for j in range(len(normalizedJukesMatrix[i])):
                    correctedMatrixFile.write(str(normalizedJukesMatrix[i][j]))
                    correctedMatrixFile.write(" ")
                correctedMatrixFile.write("\n")
            print("Matrix written to correction_output/correctedJukesMatrix.txt")
            correctedMatrixFile.close()
        elif user_model_input == "2":
            K80Matrix = correctorInstance.K80()
            normalizedK80Matrix = CorrectionToolsInstance.normalize(inputMatrix= K80Matrix)
            correctedMatrixFile = open("correction_output/correctedK80Matrix.txt", "w")
            for i in range(len(normalizedK80Matrix)):
                for j in range(len(normalizedK80Matrix[i])):
                    correctedMatrixFile.write(str(normalizedK80Matrix[i][j]))
                    correctedMatrixFile.write(" ")
                correctedMatrixFile.write("\n")
            print("Matrix written to correction_output/correctedK80Matrix.txt")
            correctedMatrixFile.close()
        elif user_model_input == "3":
            K81Matrix = correctorInstance.K81()
            normalizedK81Matrix = CorrectionToolsInstance.normalize(inputMatrix= K81Matrix)
            correctedMatrixFile = open("correction_output/correctedK81Matrix.txt", "w")
            for i in range(len(normalizedK81Matrix)):
                for j in range(len(normalizedK81Matrix[i])):
                    correctedMatrixFile.write(str(normalizedK81Matrix[i][j]))
                    correctedMatrixFile.write(" ")
                correctedMatrixFile.write("\n")
            print("Matrix written to correction_output/correctedK81Matrix.txt")
            correctedMatrixFile.close()
        elif user_model_input == "4":
            T92Matrix = correctorInstance.T92()
            normalizedT92Matrix = CorrectionToolsInstance.normalize(inputMatrix= T92Matrix)
            #Writing the corrected Matrix to an output file
            correctedMatrixFile = open("correction_output/correctedT92Matrix4.txt", "w")
            for i in range(len(normalizedT92Matrix)):
                for j in range(len(normalizedT92Matrix[i])):
                    correctedMatrixFile.write(str(normalizedT92Matrix[i][j]))
                    correctedMatrixFile.write(" ")
                correctedMatrixFile.write("\n")
            print("Matrix written to correction_output/correctedT92Matrix4.txt")
            correctedMatrixFile.close()
        elif user_model_input == "5":
            TN93Matrix = correctorInstance.TN93()
            normalizedTN93Matrix = CorrectionToolsInstance.normalize(inputMatrix= TN93Matrix)
            #Writing the corrected Matrix to an output file
            correctedMatrixFile = open("correction_output/correctedTN93Matrix.txt", "w")
            for i in range(len(normalizedTN93Matrix)):
                for j in range(len(normalizedTN93Matrix[i])):
                    correctedMatrixFile.write(str(normalizedTN93Matrix[i][j]))
                    correctedMatrixFile.write(" ")
                correctedMatrixFile.write("\n")
            print("Matrix written to correction_output/correctedTN93Matrix.txt")
            correctedMatrixFile.close()
        else:
            print("Invalid input. Please try a number 1 - 5.")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    main()