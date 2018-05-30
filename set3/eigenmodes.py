import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
import math

def laplacian_matrix(gridpoints, reduction = True):
    size = gridpoints**2
    right = np.diag(np.ones(size-1), 1)
    left = np.diag(np.ones(size-1), -1)
    up = np.diag(np.ones(size-gridpoints), gridpoints)
    down = np.diag(np.ones(size-gridpoints), -gridpoints)
    diag = -4*np.diag(np.ones(size))

    sq = math.sqrt(size)
    matrix = right + left + up + down + diag
    for i in range(size):
        for j in range(size):
            if(i <= sq or i >= size - sq-1) or i%gridpoints ==gridpoints-1 or i%gridpoints == 0:
                matrix[i,j] = 0

            elif (j <= sq or j >= size - sq-1) or j%gridpoints ==gridpoints-1 or j%gridpoints == 0:
                matrix[i,j] = 0

    if reduction:
        matrix = matrix[gridpoints+1:size-gridpoints-1, gridpoints+1:size-gridpoints-1]

        rows = []# Rows to delete
        for i in range(len(matrix[0,:])):
            if np.all(matrix[i,:] == 0):
                rows.append(i)

        while len(rows)> 0:
            row_to_remove = rows[-1]
            matrix = np.delete(matrix, row_to_remove, 0)
            matrix = np.delete(matrix, row_to_remove, 1)
            rows.remove(row_to_remove)

    return matrix


if __name__ == '__main__':
    griddimension = 40
    l = 1
    dx = l/float(griddimension)
    reduction = True

    edgesgrid = np.zeros((griddimension, griddimension))
    for i in range(griddimension):
        for j in range(griddimension):
            if i%griddimension != 0 and i%griddimension != griddimension-1 and j%griddimension != 0 and j%griddimension != griddimension -1:
                edgesgrid[i,j] = 1

    Matrix1 = laplacian_matrix(griddimension, reduction)/dx**2
    Eigenvalues, EigenVectors = linalg.eigh(Matrix1)
    print(Matrix1)
    print(Eigenvalues)
    index = 0
    for Eig in Eigenvalues:
        grid = np.zeros((griddimension,griddimension))
        k=0
        for i in range(griddimension):
            for j in range(griddimension):
                if reduction:
                    if edgesgrid[i,j] == 1:
                        grid[i,j] = EigenVectors[:,index][k]
                        k+=1
        index+=1

        if abs(Eig) < 110:
            plt.imshow(grid, origin = 'lower')
            plt.title('Frequency = ' + str(abs(round(Eig,2))))
            plt.colorbar()
            plt.show()
