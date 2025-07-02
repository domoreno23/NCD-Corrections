'''Interface for optimizing perameters of an evolutionary model
1. Uses arbitrary gusses for intial perameter values
2. Utilizes a MLE calculation to optimize perameters
3. ONLY returns optimized perameters to be used by the Corrector class
'''
import numpy as np
import math
from scipy.optimize import minimize

class PeramterOptimizer:
    def __init__(self, originalDistances, sequenceFile,):
        #These two will be passed in from the Corrector class
        self.originalDistances = originalDistances
        self.sequencesForFile = sequenceFile
        self.optimizedPerameters = []
        self.HKY85Matrix = []
        self.T92Matrix = []
    
    def HKY85(self, perams):
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
            
            k = 0.5
            
            b = 1 / (2 * (freqA + freqG) + (freqC + freqT) + 2*k * ((freqA*freqG) + (freqC * freqT)))
            
            
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
                        
                        corrList.append(x)
                self.HKY85Matrix.append(corrList)

            
            return self.HKY85Matrix
    
    def T92(self, perams):

        labels = [row[0] for row in self.originalDistances]
        data = [row[1:] for row in self.originalDistances]
        T92InitialMatrix = np.array(data, dtype=float)
        
        p_initial, q_initial, G_and_C = perams[0], perams[1], perams[2]
        
        h = 2 * G_and_C * (1 - G_and_C)
        firstLog = abs(1 - (p_initial/h) - q_initial)
        secondLog = abs(1 - (2*q_initial))
        if firstLog < 0:
            firstLog = 0.01
        if secondLog < 0:
            secondLog = 0.01
            
        d = -h * math.log(firstLog) - 1/2 * (1 - h) * math.log(secondLog)
        
        self.T92Matrix = T92InitialMatrix * d
        
        return np.array(self.T92Matrix)
    
    def TN93(self, perams):
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
        
    def GTR(self, perams):
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
    
    
    #Calculates differences between model and observed distances, core of MLE
    def loss(self, modelNumber, perams):
        #NOTE: perams are passed through here first and then passed again to each of the model functions
        model_distances = np.array([])
        labels = [row[0] for row in self.originalDistances]
        data = [row[1:] for row in self.originalDistances]
        observed_distances = np.array(data, dtype=float)
        if modelNumber == 1:
            model_distances = self.HKY85(perams)
        elif modelNumber == 2:
            model_distances = self.T92(perams)
        elif modelNumber == 3:
            model_distances = self.TN93(perams)
        elif modelNumber == 4:
            model_distances = self.GTR(perams)
        
        
        return np.sum((observed_distances - model_distances) ** 2)
    
    
    #Optimizes perameters by utilizing the differences and initial values
    def optimize(self, initialPerams, modelNumber, perameterBounds):
        initialValues = np.array(initialPerams)
        result = minimize(lambda perams: self.loss(modelNumber, perams), x0=initialValues, bounds=perameterBounds)
        
        optimizedPerams = result.x
        return optimizedPerams