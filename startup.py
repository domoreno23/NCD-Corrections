import CorrectionTools
import Corrector

def main():
    print("Welcome to the NCD corrections program!")
    
    '''Getting the user matrix file'''
    
    user_matrix_file = input("Please enter the path to the distance matrix file: ")
    try:
        with open(user_matrix_file, 'r') as file:
            # Inputting the file into an instance of CorrectionTools
            inputMatrix = CorrectionTools.CorrectionTools(user_matrix_file).createMatrix()
            
    except FileNotFoundError:
        print("Error: File not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
        
    
    '''Getting user FASTA file'''
    
    user_FASTA_file = input("Please enter the path to the FASTA file: ")
    try:
        with open(user_FASTA_file, 'r') as file:
            correctionInstance = Corrector(inputMatrix, user_FASTA_file)
    except FileNotFoundError:
        print("Error: File not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    #Creating the corrected distance matrix that the user had input
    print("Please choose a correction model:")
    user_model_input = input("Jukes Cantor (1), F81 (2), HKY85 (3)")
    
    
    try:
        if user_model_input == "1":
            jukesMatrix = correctionInstance.jukesCantor()
            print(jukesMatrix)
        elif user_model_input == "2":
            F81Matrix = correctionInstance.F81()
            print(F81Matrix)
        elif user_model_input == "3":
            HKy85Matrix = correctionInstance.HKY85()
            print(HKy85Matrix)
        else:
            print("Invalid input. Please try a number 1 - 3.")
    except Exception as e:
        print(f"An error occurred: {e}")


'''
if __name__ == '__main__':

    tools = CorrectionTools.CorrectionTools('LargeListMatrixMfc2_37_0_0.txt')
    tools.createMatrix()
    dists = tools.distanceMatrix
    print(dists)

    corr = Corrector.Corrector(dists, "gene447_final_unaligned.txt")

    ### Below are the calls to the correction methods
    ###Jukes Cantor and F81 work. I can't remember if I ever got F81Mod running

    # corr.jukesCantor()
    # jukes = corr.JukesCantorMatrix
    # # print(jukes)
    #
    # corr.F81()
    # f81 = corr.F81Matrix
    # print(f81)
    # print(corr.frequenciesPerSequence)
    # print(len(corr.frequenciesPerSequence))

    # corr.F81Mod()
    # f81Mod = corr.F81ModMatrix
    # print(f81Mod)
'''