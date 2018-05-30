import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
import math

def laplacian_matrix(width, length, reduction = True):
    size = width*length
    right = np.diag(np.ones(size-1), 1)
    left = np.diag(np.ones(size-1), -1)
    up = np.diag(np.ones(size-width), width)
    down = np.diag(np.ones(size-width), -width)
    diag = -4*np.diag(np.ones(size))

    matrix = right + left + up + down + diag
    for i in range(size):
        for j in range(size):
            if(i <= width or i >= size - width-1) or i%width ==width-1 or i%width == 0:
                matrix[i,j] = 0

            elif (j <= width or j >= size - width-1) or j% length == length-1 or j%length == 0:
                matrix[i,j] = 0

    if reduction:
        rows = []# Rows to delete
        collumns = []
        for i in range(len(matrix[0,:])):
            if np.all(matrix[i,:] == 0):
                rows.append(i)

        while len(rows)> 0:
            row_to_remove = rows[-1]
            matrix = np.delete(matrix, row_to_remove, 0)
            rows.remove(row_to_remove)

        for i in range(len(matrix[0, :])):
            if sum(matrix[:, i]) == 1 or np.all(matrix[:, i] == 0):
                collumns.append(i)

        while len(collumns) > 0:
            column_to_remove = collumns[-1]
            matrix = np.delete(matrix, column_to_remove, 1)
            collumns.remove(column_to_remove)
            
    print(matrix)
    return matrix


if __name__ == '__main__':
    width = 40
    length = 80
    l=1
    dx = l/float(width)
    reduction = True

    edgesgrid = np.zeros((width, length))
    for i in range(width):
        for j in range(length):
            if i % width != 0 and i % width != width -1 and j % length != 0 and j % length != length -1:
                edgesgrid[i,j] = 1

    Matrix1 = laplacian_matrix(width, length, reduction)/dx**2
    Eigenvalues, EigenVectors = linalg.eig(Matrix1)

    index = 0
    for Eig in Eigenvalues:
        grid = np.zeros((width,length))
        k=0
        for i in range(width):
            for j in range(length):
                if reduction:
                    if edgesgrid[i,j] == 1:
                        grid[i,j] = EigenVectors[:,index][k]
                        k+=1
        index+=1
        
        if abs(Eig) <120:
            plt.imshow(grid, origin = 'lower')
            plt.title('Frequency = ' + str(abs(round(Eig,2))))
            plt.colorbar()
            plt.show()
