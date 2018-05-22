import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
import math

def laplacian_matrix(size):
    # L = 1
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

    index = 0
    while index < size:
        if(index % sq == 0 or index % sq == 1):
            for j in range(size):
                matrix[index, j] = 0
        index += 1

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
    n = 5
    m = laplacian_matrix(n**2)
    eigv = linalg.eig(m)
    minv = np.argmin(eigv[0])

    xticks = ['v'+str(i)+str(j) for i in xrange(1,n+1) for j in xrange(1,n+1)]
    x = np.arange(len(eigv[1][minv]))
    print x
    plt.figure()
    plt.bar(x, eigv[1][minv])
    plt.xticks(x,xticks, rotation=90)
    plt.title("Eigenvector values, frequency = " + str(eigv[0][minv]))
    plt.show()

    # print (circle_matrix(16))
