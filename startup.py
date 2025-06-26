import CorrectionTools
import Corrector

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
        print("Error: File not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    #Creating the corrected distance matrix that the user had input
    print("Please choose a correction model:")
    user_model_input = input("Jukes Cantor (1), F81 (2), HKY85 (3)")
    
    
    #Takes instance of the inputted matrix and runs the selected correction
    #NOTE: Implement feature that allows to have multiple corrections within a session
    try:
        if user_model_input == "1":
            
            jukesMatrix = correctorInstance.jukesCantor()
            normalizedJukesMatrix = CorrectionToolsInstance.normalize(inputMatrix= jukesMatrix)
            
            #Writing the corrected Matrix to an output file
            correctedMatrixFile = open("correctedJukesMatrix.txt", "w")
            
            for i in range(len(normalizedJukesMatrix)):
                for j in range(len(normalizedJukesMatrix[i])):
                    correctedMatrixFile.write(str(normalizedJukesMatrix[i][j]))
                    correctedMatrixFile.write(" ")
                correctedMatrixFile.write("\n")
            print("Matrix written to correctedJukesMatrix.txt")
            correctedMatrixFile.close()
        elif user_model_input == "2":
            F81Matrix = correctorInstance.F81()
            print(F81Matrix)
        elif user_model_input == "3":
            HKy85Matrix = correctorInstance.HKY85()
            print(HKy85Matrix)
        else:
            print("Invalid input. Please try a number 1 - 3.")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    main()