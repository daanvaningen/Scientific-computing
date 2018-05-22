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
    # for i in range(size):
    #     for j in range(size):
    #         if(i <= sq or i >= size - sq-1):
    #             matrix[i,j] = 0
    #         elif(j <= sq or j >= size - sq-1):
    #             matrix[i,j] = 0

    return matrix

def pythagoras(i, j, radius):
    return math.sqrt((radius - float(i))**2 + (radius - float(j))**2)

def circle_matrix(L):
    matrix = np.zeros((L, L))
    radius = 0.5*L
    for i in range(L):
        skip = False
        if pythagoras(i,i, radius) > radius:
            skip = True
        places = [i-4, i-1, i+1, i+4]
        for j in range(L):
            if not skip:
                if(i == 0 or i == L-1 or j == 0 or j == L-1):
                    matrix[i,j] = 0
                elif(pythagoras(i,j,radius) > radius):
                    matrix[i,j] = 0
                elif j in places:
                    matrix[i,j] = 1
                elif j == i:
                    matrix[i,j] = -4

    return matrix

if __name__ == '__main__':
    # print laplacian_matrix(16)
    print (circle_matrix(16))
