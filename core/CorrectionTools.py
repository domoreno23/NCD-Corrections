'''Interface for creating and normalizing a distance matrix'''
import numpy as np
from sklearn.preprocessing import MinMaxScaler

class CorrectionTools:

    def __init__(self, distanceMatrixFile):
        distMatrixFile = open(distanceMatrixFile)
        self.rawDists = distMatrixFile.readlines()
        self.distanceMatrix = None
        self.normalizedMatrix = None
        
        
    
    def createMatrix(self): #preprocessing
        newMatrix = []
        for x in range(1, len(self.rawDists)):
            dsts = self.rawDists[x].split()
            # Keep the label as is, convert the rest to float
            row = [dsts[0]] + [float(val) for val in dsts[1:]]
            newMatrix.append(row)
        self.distanceMatrix = newMatrix
        return self.distanceMatrix
    
    
    
    def normalize(self, inputMatrix): #preprocessing
        # Extract only the numeric part for normalization
        labels = [row[0] for row in inputMatrix]
        data = [row[1:] for row in inputMatrix]
        
        normalizer = MinMaxScaler()
        normalized_data = normalizer.fit_transform(data)
        normalizedMatrix = [[label] + list(row) for label, row in zip(labels, normalized_data)]
        self.normalizedMatrix = normalizedMatrix
        return self.normalizedMatrix