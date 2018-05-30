import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig
import math

def laplacian_matrix(width, length, Vector, reduction = True):
    size = width*length
    right = np.diag(np.ones(size-1), 1)
    left = np.diag(np.ones(size-1), -1)
    up = np.diag(np.ones(size-width), width)
    down = np.diag(np.ones(size-width), -width)
    diag = -4*np.diag(np.ones(size))

    V_vector = Vector
    rows_to_del = []
    matrix = right + left + up + down + diag

    if reduction:
        k=0
        for i in V_vector:
            if i == 0:
                rows_to_del.append(k)
                
            k+=1
            
        while len(rows_to_del) > 0:
            last_element = rows_to_del[-1]
            matrix = np.delete(matrix, last_element, 0)
            matrix = np.delete(matrix, last_element, 1)

            rows_to_del.remove(last_element)
            
    print(matrix)
    return matrix


if __name__ == '__main__':
    width = 40             # Aantal grid points in breedte
    length = 40             # Als dit een cirle is, dan naar 40 zetten.
    l=1                     # lengte van drum, als je een circle van straal 1 wilt moet je hier 2 invoeren.
    dx = l/float(width)
    reduction = True        # Whether or not you want the matrix to be reducted, all zeros are gone
    circle = False          # Whether the domain is a circle or not
    time_dependent = True
    
    if not circle:
        edgesgrid = np.zeros((width, length))
        for i in range(width):
            for j in range(length):
                if i % width != 0 and i % width != width -1 and j % length != 0 and j % length != length -1:
                    edgesgrid[i,j] = 1
    else:
        edgesgrid = np.zeros((width, width))
        for i in range(width):
            for j in range(width):
                if math.sqrt((i*dx*l-l)**2 + (j*dx*l-l)**2) < l:
                    edgesgrid[i,j] = 1
    edgesgrid = edgesgrid.flatten()

    Matrix1 = laplacian_matrix(width, length, edgesgrid)/dx**2
    Eigenvalues, EigenVectors = eig(Matrix1)

    index = 0
    for Eig in Eigenvalues:
        grid = np.zeros((width,length))
        k=0
        m=0
        for i in range(width):
            for j in range(length):
                if reduction:
                    if edgesgrid[m] == 1:
                        grid[i,j] = EigenVectors[:,index][k]
                        k+=1
                    m+=1
        index+=1
        
        if abs(Eig) <120:
            if time_dependent:
                for t in range(0,10):
                    grid = grid * (math.sin(abs(Eig)*t/10.0) + math.cos(abs(Eig)*t/10.0))

                    plt.imshow(grid, origin = 'lower')
                    plt.title('$\lambda$ = ' + str(abs(round(Eig,2))) + '  and t = ' + str(t/10.0))
                    plt.colorbar()
                    plt.show()
