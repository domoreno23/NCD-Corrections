# this class provides necessary tools for pre- and post- data processing

class CorrectionTools:

#Not being used yet; was in the process of building this out

    def __init__(self, distanceMatrixFile):
        distMatrixFile = open(distanceMatrixFile)
        self.rawDists = distMatrixFile.readlines()
        self.distanceMatrix = None
        # self.normalizedMatrix = self.normalize(self.distanceMatrix)


    def createMatrix(self): #preprocessing
        newMatrix = []
        newMatrix.append(self.rawDists[0].strip())
        for x in range(1, len(self.rawDists)):
            dsts = self.rawDists[x].split()
            newMatrix.append(dsts)
        self.distanceMatrix = newMatrix
        return self.distanceMatrix

    # def normalize(self): #preprocessing
    #     # need to determine if this is necessary for all corrections or just jukes
    #     # is there a better library I can use for this?? unsure why I wrote my own code for this?
    #
    #     maxDst = max([max(matrix[i][1:]) for i in range(1, len(matrix))])
    #     # print(maxDst)
    #     minDst = 0.0
    #     newMatrix = []
    #     newMatrix.append(matrix[0])
    #     for i in range(1, len(matrix)):  # iterate thru rows
    #         row = []
    #         row.append(matrix[i][0])  # species name
    #         for j in range(1, len(matrix[i])):  # iterate thru columns (dists)
    #             dist = float(matrix[i][j])
    #             normDist = (((dist - minDst) / (float(maxDst) - minDst)))
    #             row.append(normDist)
    #         newMatrix.append(row)
    #     return newMatrix
