'''Calculates the Euclidean distance between all points between two matrices.'''
import numpy as np
import core.CorrectionTools as CorrectionTools

def euclidean_distance(matrix1, matrix2):
  CorrectionToolsInstanceOne = CorrectionTools.CorrectionTools(matrix1)
  matrix1 = CorrectionToolsInstanceOne.createMatrix()
  CorrectionToolsInstanceTwo = CorrectionTools.CorrectionTools(matrix2)
  matrix2 = CorrectionToolsInstanceTwo.createMatrix()
  
  labels_1 = [row[0] for row in matrix1]
  data_1 = [row[1:] for row in matrix1]
  
  labels_2 = [row[0] for row in matrix2]
  data_2 = [row[1:] for row in matrix2]
  
  data_1 = np.array(data_1, dtype=float)
  data_2 = np.array(data_2, dtype=float)
  
  if data_1.shape != data_2.shape:
    print(data_1.shape, data_2.shape)
    raise ValueError("Matrices must have the same dimensions for Euclidean distance calculation.")
  
  distance = np.sqrt(np.sum((data_1 - data_2) ** 2))/(data_1.shape[0] ** 2 - data_1.shape[0])
  
  
  return distance

distance = euclidean_distance(matrix1="test_files/LargeListMatrix7Zip_153_0_0.txt", matrix2="correction_output/correctedJukesMatrix.txt")

print(f"The Euclidean distance between the two matrices is: {distance}")