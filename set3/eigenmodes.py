import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

class Grid:
    def __init__(self, size):
        self.size = size
        self.grid = np.zeros((size+2, size+2))


def laplacian_matrix(size):
    up = -np.diag(np.ones(size-1), 1)
    low = -np.diag(np.ones(size-1), -1)
    diag = 2*np.diag(np.ones(size))

    return up + low + diag

if __name__ == '__main__':
    grid = Grid(10)

    print laplacian_matrix(10)
