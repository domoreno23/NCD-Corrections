'''Interface for optimizing perameters of an evolutionary model
1. Uses arbitrary gusses for intial perameter values
2. Utilizes a MLE calculation to optimize perameters
3. ONLY returns optimized perameters to be used by the Corrector class
'''
import numpy as np
import math
from scipy.optimize import minimize

class PeramterOptimizer:
  def __init__(self, originalDistances, sequenceFile, frequenciesForFile):
    #These two will be passed in from the Corrector class
    self.originalDistances = originalDistances
    self.sequenceFile = sequenceFile
    self.frequenciesForFile = frequenciesForFile
    self.optimizedPerameters = []
    self.HKY85Matrix = []
    self.T92Matrix = []
  
  def HKY85(self, perams, *args):
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
  
  def T92(self, perams, *args):
    self.T92Matrix.append(self.originalDistances[0])
    ##Setting perameters for T92
    freqG = self.frequenciesForFile["G"] / self.frequenciesForFile["Total"]
    freqC = self.frequenciesForFile["C"] / self.frequenciesForFile["Total"]
    
    
    ##NOTE: normalize frequencies before production use
    G_and_C = (freqG + freqC)/2
    
    #where p is the proportion of sites that show transitional differences and q is the proportion of sites that show transversional differences.
    p_guess = 0.1
    q_guess = 0.1
    h = 2 * G_and_C * (1 - G_and_C)
    d = -h * math.log(1 - (p_guess/h) - q_guess) - 1/2 * (1 - h) * math.log(1-2*q_guess)
    
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
                corrList.append(d)
        self.T92Matrix.append(corrList)
    
    return self.T92Matrix
  
  def TN93(self, perams, *args):
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
    
  def GTR(self, perams, *args):
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
  def loss(self, modelNumber, perams, observed_distances):
    #NOTE: None of the models have perameters because all of their initial values are arbitrarily initialized within their respective functions
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
  def optimize(self, initialPerams, modelNumber):
    initialValues = np.array(initialPerams)
    bounds = (0, 1)
    
    result = minimize(self.loss(modelNumber, initialValues, self.originalDistances), x0=initialValues, bounds=bounds)
    
    optimizedPerams = result.x
    return optimizedPerams