# Goals (class Corrector):
# - implements the called corrections on a given distance matrix that has been preprocessed through Correction Tools
# - class variable for original distanceMatrix
# - class variable to store each corrected matrix (one for each correction)
# - functions for each correction

import math
import numpy as np 
class Corrector:

    def __init__(self, distanceMatrix, sequenceFile = None):
        self.originalDistances = distanceMatrix
        self.JukesCantorMatrix = []
        self.F81Matrix = []
        self.HKY85Matrix = []
        self.T92Matrix = []
        self.TN93Matrix = []
        self.GTRMatrix = []


        #Needed for F81 and F81 MOD
        #Counts frequencies for each of the nucleotides if the file exists
        if sequenceFile is not None:
            self.F81Matrix = []
            self.F81ModMatrix = []
            seqFile = open(sequenceFile)
            content2 = seqFile.readlines()
            self.sequences = content2
            self.frequenciesForFile = {"A" : 0, "T": 0, "G": 0, "C": 0, "Total": 0}
            self.frequenciesPerSequence = [] #list of dicts: each inner dict i contains the ATGC freqs for the associated sequence i
            self.findFrequencies()
            
            
    ##Will be used for methods that require the use of transversion and transition rates
    def __calculateP_and_Q(self):
        
        pass


    def findFrequencies(self):
        # We need total nucleotides and totalA,T,G,C for the entire sequence file for F81
        # We need total nucleos and total A,T,G,C for each individual sequence as well for F81 Mod
        # ASK: what is N in the unaligned sequence? do we include that in the seq Length??

        for s in range(1, len(self.sequences), 2):
            #seqLength = len(self.sequences[s])
            #keeps track of totals for each individual sequence (could be put in dictionary to be cleaner)
            freqsForSequence = {"A" : 0, "T": 0, "G": 0, "C": 0, "Total": 0}
            # totalNucleo = 0
            # freqA = 0
            # freqT = 0
            # freqG = 0
            # freqC = 0
            for ch in self.sequences[s]:
                if ch in ("A","T","G","C"):
                    self.frequenciesForFile[ch] += 1
                    self.frequenciesForFile["Total"] += 1
                    freqsForSequence[ch] += 1
                    freqsForSequence["Total"] += 1
            self.frequenciesPerSequence.append(freqsForSequence)
            
    

    def jukesCantor(self):
    # Separate labels and numeric data
        labels = [row[0] for row in self.originalDistances]
        data = [row[1:] for row in self.originalDistances]
        jukesOriginal = np.array(data, dtype=float)
        # Apply the Jukes-Cantor correction
        corrected = (-0.75 * np.log(1 - (jukesOriginal * (3/4))))
        # Combine labels back with corrected data
        self.JukesCantorMatrix = [[label] + list(row) for label, row in zip(labels, corrected)]
        return self.JukesCantorMatrix
        
        '''
    def jukesCantor(self):
        ##print(self.originalDistances)
        # corrMatrix = []
        # print(normDistMatrix)
        self.JukesCantorMatrix.append(self.originalDistances[0])
        ##
        # normalize this matrix based on its min/max ---  between 0 and 1
        ##
        # print(normDistMatrix)
        for i in range(1, len(self.originalDistances)):
            # print(normDistMatrix[i])
            corrList = []
            corrList.append(self.originalDistances[i][0])
            for j in range(1, len(self.originalDistances[i])):
                dist = float(self.originalDistances[i][j])
                if dist == 0.0:
                    # normDistMatrix[i][j] = dist
                    corrList.append(dist)
                else:  # correction
                    ##
                    # Since distances are normalized
                    # ==> multiply by .74 to ensure range [0-.74] for JUKES_CANTOR requirement
                    ##
                    x = 1 - (dist * .74 / .75)
                    # print(dist)
                    # normDistMatrix[i][j] = (-.75 * (math.log(x)))
                    corrList.append((-.75 * (math.log(x))))
            self.JukesCantorMatrix.append(corrList)
        
        return self.JukesCantorMatrix
        '''
    def F81(self):

        #Need frequencies for whole file
        freqA = self.frequenciesForFile["A"] / self.frequenciesForFile["Total"]
        freqT = self.frequenciesForFile["T"] / self.frequenciesForFile["Total"]
        freqG = self.frequenciesForFile["G"] / self.frequenciesForFile["Total"]
        freqC = self.frequenciesForFile["C"] / self.frequenciesForFile["Total"]

        b = 1 - (freqA ** 2) - (freqT ** 2) - (freqG ** 2) - (freqC ** 2)
        print(b)

        # apply correction
        self.F81Matrix.append(self.originalDistances[0])
        ##
        # normalize this matrix based on its min/max ---  between 0 and 1
        ##
        # print(normDistMatrix)
        for i in range(1, len(self.originalDistances)):
            # print(normDistMatrix[i])
            corrList = []
            corrList.append(self.originalDistances[i][0])
            for j in range(1, len(self.originalDistances[i])): # j indicies for dists go from [1 to 37]
                dist = float(self.originalDistances[i][j])
                if dist == 0.0:
                    corrList.append(dist)
                else:  # correction
                    ##
                    # Since distances are normalized
                    # ==> multiply by .74 to ensure range [0-.74] for JUKES_CANTOR requirement
                    ##
                    x = 1 - (dist * .74 / b)  # maybe dont want to normalize to .74 for future use since b might not be .75 (for now prob fine since b is abt .74)
                    # print(x)
                    # print(dist)
                    # normDistMatrix[i][j] = (-.75 * (math.log(x)))
                    corrList.append((-b * (math.log(x))))
            self.F81Matrix.append(corrList)

        # normCorrDist = normalize(corrMatrix)
        # return normCorrDist
        return self.F81Matrix


    def F81Mod(self):

        # find min B then normalize ??? or just do .74
        self.F81ModMatrix.append(self.originalDistances[0])
        for i in range(1, len(self.originalDistances)):  # correction
            correctedRow = []
            correctedRow.append(self.originalDistances[i][0])
            # need freq of atgc for i
            freqListI = self.frequenciesPerSequence[i-1]  # {A:, T:, G:, C:, Total:}
            freqAI = freqListI["A"] / freqListI["Total"]
            freqTI = freqListI["T"] / freqListI["Total"]
            freqGI = freqListI["G"] / freqListI["Total"]
            freqCI = freqListI["C"] / freqListI["Total"]
            for j in range(1, len(self.originalDistances[i])):
                dist = float(self.originalDistances[i][j])
                if dist == 0.0:
                    correctedRow.append(dist)
                else:
                    # get b1:
                    freqA = self.frequenciesForFile["A"] / self.frequenciesForFile["Total"]
                    freqT = self.frequenciesForFile["T"] / self.frequenciesForFile["Total"]
                    freqG = self.frequenciesForFile["G"] / self.frequenciesForFile["Total"]
                    freqC = self.frequenciesForFile["C"] / self.frequenciesForFile["Total"]
                    b1 = 1 - (freqA ** 2) - (freqT ** 2) - (freqG ** 2) - (freqC ** 2)

                    # get b2:
                    freqListJ = self.frequenciesPerSequence[j-1]
                    # {A:, T:, G:, C:, Total:} #j-1 bc this list goes from 0 to 36
                    # and the dist indices (j) are from 1 to 37 so j corresponds to j-1 in the freqList
                    freqAJ = freqListJ["A"] / freqListJ["Total"]
                    freqTJ = freqListJ["T"] / freqListJ["Total"]
                    freqGJ = freqListJ["G"] / freqListJ["Total"]
                    freqCJ = freqListJ["C"] / freqListJ["Total"]

                    b2 = 1 - ((freqAI * freqAJ) ** 2) - ((freqTI * freqTJ) ** 2) - ((freqGI * freqGJ) ** 2) - (
                            (freqCI * freqCJ) ** 2)

                    x = 1 - ((dist * .74) / b2)
                    # print((-b1 * (math.log(x))))
                    correctedRow.append((-b1 * (math.log(x))))
            self.F81ModMatrix.append(correctedRow)

        # print("f81")
        # print(newMatrix)
        # normMatrix = normalize(newMatrix)
        # print(normMatrix)
        # return (normMatrix)
        return self.F81ModMatrix
    
    '''HKY85, distinguishes between the rate of transitions and transversions (using the κ parameter), and it allows unequal base frequencies'''
    def HKY85(self):
        #Show user the original distance matrix
        print(self.originalDistances)
        
        # corrMatrix = []
        # print(normDistMatrix)
        self.HKY85Matrix.append(self.originalDistances[0])
        ##
        # normalize this matrix based on its min/max ---  between 0 and 1
        ##
        # print(normDistMatrix)
        #Need frequencies for whole file
        freqA = self.frequenciesForFile["A"] / self.frequenciesForFile["Total"]
        freqT = self.frequenciesForFile["T"] / self.frequenciesForFile["Total"]
        freqG = self.frequenciesForFile["G"] / self.frequenciesForFile["Total"]
        freqC = self.frequenciesForFile["C"] / self.frequenciesForFile["Total"]
        
        k = 0
        
        b = 1 / (2 * (freqA + freqG) + (freqC + freqT) + 2*k * ((freqA*freqG) + (freqC * freqT)))
        print(b)
        
        for i in range(1, len(self.originalDistances)):
            # print(normDistMatrix[i])
            corrList = []
            corrList.append(self.originalDistances[i][0])
            for j in range(1, len(self.originalDistances[i])):
                dist = float(self.originalDistances[i][j])
                if dist == 0.0:
                    # normDistMatrix[i][j] = dist
                    corrList.append(dist)
                else:  # correction
                    ##
                    # Since distances are normalized
                    ##
                    x = (dist * b)
                    # print(dist)
                    # normDistMatrix[i][j] = (-.75 * (math.log(x)))
                    corrList.append(x)
            self.HKY85Matrix.append(corrList)

        ### Re-normalize correct distances to [0-1]
        # print("HKY85 Corrections")
        # print(corrMatrix)
        # normCorrDist = normalize(corrMatrix)
        # print(normCorrDist)
        # return normCorrDist
        return self.HKY85Matrix
    
    
    ##NOTE: ask if we could multiply matricies in a while loop until the multiplaction converges??
    def T92(self):
        print(self.originalDistances)
        
        self.T92Matrix.append(self.originalDistances[0])
        ##Setting perameters for T92
        freqG = self.frequenciesForFile["G"] / self.frequenciesForFile["Total"]
        freqC = self.frequenciesForFile["C"] / self.frequenciesForFile["Total"]
        
        
        ##NOTE: normalize frequencies before production use
        G_and_C = (freqG + freqC)/2
        
        #where p is the proportion of sites that show transitional differences and q is the proportion of sites that show transversional differences.
        p = 0
        q = 0
        h = 2 * G_and_C * (1 - G_and_C)
        d = -h * math.log(1 - (p/h) - q) - 1/2 * (1 - h) * math.log(1-2*q)
        
        for i in range(1, len(self.originalDistances)):
            # print(normDistMatrix[i])
            corrList = []
            corrList.append(self.originalDistances[i][0])
            for j in range(1, len(self.originalDistances[i])):
                dist = float(self.originalDistances[i][j])
                if dist == 0.0:
                    # normDistMatrix[i][j] = dist
                    corrList.append(dist)
                else:  # correction
                    ##
                    # Since distances are normalized
                    ##
                    x = 0
                    # print(dist)
                    # normDistMatrix[i][j] = (-.75 * (math.log(x)))
                    corrList.append(x)
            self.T92Matrix.append(corrList)
            
        return self.T92Matrix
    
    def TN93(self):
        print(self.originalDistances)
        
        self.TN93Matrix.append(self.originalDistances[0])
        ##Setting perameters for TN93
        #Purines rate (A and G)
        freqA = self.frequenciesForFile["A"] / self.frequenciesForFile["Total"]
        freqG = self.frequenciesForFile["G"] / self.frequenciesForFile["Total"]
        purine_transition_rate = freqA + freqG / self.frequenciesForFile["Total"]
        
        #Pyrimidines rate (C and T)
        freqC = self.frequenciesForFile["C"] / self.frequenciesForFile["Total"]
        freqT = self.frequenciesForFile["T"] / self.frequenciesForFile["Total"]
        pyrimidine_transition_rate = freqC + freqT / self.frequenciesForFile["Total"]
        
        for i in range(1, len(self.originalDistances)):
            # print(normDistMatrix[i])
            corrList = []
            corrList.append(self.originalDistances[i][0])
            for j in range(1, len(self.originalDistances[i])):
                dist = float(self.originalDistances[i][j])
                if dist == 0.0:
                    # normDistMatrix[i][j] = dist
                    corrList.append(dist)
                else:  # correction
                    ##
                    # Since distances are normalized
                    ##
                    x = 0
                    # print(dist)
                    # normDistMatrix[i][j] = (-.75 * (math.log(x)))
                    corrList.append(x)
            self.TN93Matrix.append(corrList)
            
        return self.TN93Matrix
    
    def GTR(self):
        print(self.originalDistances)
        
        self.GTRMatrix.append(self.originalDistances[0])
        ##Setting perameters for GTR
        
        ##Rate of A <-> C
        alpha = 0
        
        ## Rate of A <-> C
        beta = 0
        
        ## Rate of A <-> T 
        Gamma = 0
        
        ## Rate of G <-> C rate
        delta = 0
        
        ## Rate of G <-> T
        epsilon = 0
        
        ##Rate of C <-> T
        eta = 0
        
        for i in range(1, len(self.originalDistances)):
            # print(normDistMatrix[i])
            corrList = []
            corrList.append(self.originalDistances[i][0])
            for j in range(1, len(self.originalDistances[i])):
                dist = float(self.originalDistances[i][j])
                if dist == 0.0:
                    # normDistMatrix[i][j] = dist
                    corrList.append(dist)
                else:  # correction
                    ##
                    # Since distances are normalized
                    ##
                    x = 0
                    # print(dist)
                    # normDistMatrix[i][j] = (-.75 * (math.log(x)))
                    corrList.append(x)
            self.GTRMatrix.append(corrList)
            
        return self.GTRMatrix
    