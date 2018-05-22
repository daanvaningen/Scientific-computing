import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
import math

def laplacian_matrix(size):
    right = np.diag(np.ones(size-1), 1)
    left = np.diag(np.ones(size-1), -1)
    up = np.diag(np.ones(size-4), 4)
    down = np.diag(np.ones(size-4), -4)
    diag = -4*np.diag(np.ones(size))

    sq = math.sqrt(size)
    matrix = right + left + up + down + diag
    for i in range(size):
        for j in range(size):
            if(i <= sq or i >= size - sq-1):
                matrix[i,j] = 0
            elif(j <= sq or j >= size - sq-1):
                matrix[i,j] = 0

    return matrix

if __name__ == '__main__':
    print laplacian_matrix(16)
