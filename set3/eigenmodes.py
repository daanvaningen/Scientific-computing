import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

class Grid:
    def __init__(self, size):
        self.size = size
        self.grid = np.zeros((size+2, size+2))

def SetMatrix(n):
    # n is the dimension of grid, therefore the matrix will have dimension n^2
    M = np.zeros((n**2,n**2))
