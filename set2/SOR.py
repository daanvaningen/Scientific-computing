import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
import math

class Grid:
    def __init__(self, N, Epsilon, Max_Iterations, n, w=1, Object= 0, startingrid = False):
        self.MaxIt = Max_Iterations         # stop condition
        self.N = N                          # Grid size
        self.n = n                          # The nebula
        self.Iterations = 0                 # No. of current iteration
        self.Epsilon = Epsilon              # Epsilon (when do we accept change to be neglected.)
        if startingrid:
            self.grid = Anal_grid
        else:
            self.grid = np.zeros((N, N))        # Making grid of n x n
        self.w = w                         # The convergence measures
        self.Object = Object                # A matrix containing the information of objects
        self.const = w/4.0                  # Value dependent on W
    

    # The grid will update untill certain criteria are met.
    def run(self):
        while self.Iterations < self.MaxIt:
            self.Iterations += 1            # Keep track of iterations
            Pgrid = self.creation_Pgrid()        # Creates the grid with chances of turning
            tempgrid = np.zeros((self.N, self.N))
            for i in range(self.N):
                for j in range(self.N):
                    value = self.get_boundry_values(i, j, tempgrid, Pgrid)
                    tempgrid[i,j] = value
                  
            self.grid = tempgrid

        return self.grid, self.Object, Pgrid

    # a function that creates the Pgrid
    def creation_Pgrid(self):
        Pgrid = np.zeros((N,N))
        sum_of_concentration = 0
        for i in range(N):
            for j in range(N):
                if self.Object[i,j] == 0 and self.check_object_boundary(i,j):
                    sum_of_concentration += self.grid[i,j]**self.n
                    Pgrid[i,j] = self.grid[i,j]**self.n
        if sum_of_concentration != 0:
            Pgrid = Pgrid/float(sum_of_concentration)

        return Pgrid

    def get_boundry_values(self, i, j, tempgrid, Pgrid):
        # Checking whether it is part of an object
        if self.Object[i,j] ==1:
            return 0

        elif Pgrid[i,j] > 0:
            RandomNumber = random.random()
            if RandomNumber < Pgrid[i,j]:
                self.Object[i,j] = 1
                
                return 0
            else:
                if(j == 0 or j == self.N-1):
                    return j/(self.N-1)
                elif(i == 0):
                    return  self.const*(self.grid[1, j] +
                        self.grid[N-1, j] + self.grid[0, j+1] + tempgrid[0, j-1]) + (1-self.w)*self.grid[i,j]
                elif(i == self.N-1):
                    return self.const*(self.grid[1, j] +
                        tempgrid[i-1, j] + self.grid[i, j+1] + tempgrid[i, j-1]) + (1-self.w)*self.grid[i,j]
                else:
                    return self.const*(self.grid[i+1, j] +
                            tempgrid[i-1, j] + self.grid[i, j+1] + tempgrid[i, j-1]) + (1-self.w)*self.grid[i,j]
        else:

            if(j == 0 or j == self.N-1):
                return j/(self.N-1)
            elif(i == 0):
                return  self.const*(self.grid[1, j] +
                    self.grid[N-1, j] + self.grid[0, j+1] + tempgrid[0, j-1]) + (1-self.w)*self.grid[i,j]
            elif(i == self.N-1):
                return self.const*(self.grid[1, j] +
                    tempgrid[i-1, j] + self.grid[i, j+1] + tempgrid[i, j-1]) + (1-self.w)*self.grid[i,j]
            else:
                return self.const*(self.grid[i+1, j] +
                        tempgrid[i-1, j] + self.grid[i, j+1] + tempgrid[i, j-1]) + (1-self.w)*self.grid[i,j]

    # This function should return true of it is next to the object
    def check_object_boundary(self, i,j):
        if (i > 0 and i < N-1) and (j>0 and j < N-1):
            if self.Object[i-1,j] == 1 or self.Object[i+1, j] == 1 or self.Object[i,j-1] ==1 or self.Object[i,j+1] == 1:
                return True
            else:
                return False

        elif i == 0 and j == 0:
            if self.Object[1,0] == 1 or self.Object[0,1] == 1 or self.Object[N-2, 0]:
                return True
            else:
                return False
        elif i == 0 and j == N-1:
            if self.Object[1,N-1] == 1 or self.Object[0,N-2]==1 or self.Object[N-2, N-1] == 1:
                return True
            else:
                return False
        elif i == N-1 and j == 0:
            if self.Object[N-1,1] == 1 or self.Object[1, 0] == 1 or self.Object[N-2,0] == 1:
                return True
            else:
                return False
        elif i == N-1 and j == N-1:
            if self.Object[1, N-1] == 1 or self.Object[N-2, N-1] == 1 or self.Object[N-1, N-2] == 1:
                return True
            else:
                return False
        elif i == 0:
            if self.Object[1, j] == 1 or self.Object[0,j-1] == 1 or self.Object[0,j+1] == 1 or self.Object[N-2,j] == 1:
                return True
            else:
                return False
        elif j == 0:
            if self.Object[i+1, j] == 1 or self.Object[i-1,j] == 1 or self.Object[i, 1] == 1:
                return True
            else:
                return False

        elif i == N-1:
            if self.Object[i-1,j] == 1 or self.Object[1, j] == 1 or self.Object[i, j+1] == 1 or self.Object[i, j-1] == 1:
                return True
            else:
                return False

        elif j == N-1:
            if self.Object[i-1, j] == 1 or self.Object[i+1,j] == 1 or self.Object[i, j-1] == 1:
                return True
            else:
                return False

def ObjectCreation(size):
    grid = np.zeros((N, N))
    middle = int(N/2.0)

    grid[middle, 0] = 1

    return grid

# This function creates the grid analytically, so what it should be at the end of all iterations
# The reason why we use this is so we can start with this grid and then let the object grow into it.
def AnalyticalGridCreation(size):
    grid = np.zeros((size,size))
    Dc = 1/float(size)
    for i in range(size):
        for j in range(size):
            grid[i,j] = j * Dc

    return grid

if __name__ == '__main__':
    # Initial conditions
    N = 512                   # Grid size
    n = 2.0                    # Nebula
    w = 1.75                     # For the SOR model
    Epsilon = 0                  # The maximum change you want  to stop iterating
    Max_Iterations = 2000      # The maximum iterations before shutdown
    Line_start = True           # True means that you start with the analytical values

    # The analytical grid is the starting grid.
    Anal_grid = AnalyticalGridCreation(N)
    grid_ = ObjectCreation(N)
    SOR = Grid(N, Epsilon, Max_Iterations, n, w, grid_, Line_start)
    Equilibrium = SOR.run()

    for i in range(N):
        for j in range(N):
            if Equilibrium[1][i,j] == 1:
                Equilibrium[0][i,j] = 1.1

    plt.imshow(np.rot90(Equilibrium[0],3), origin = 'lower', cmap = "hot")
    plt.title("SOR with N=" + str(N) + ", w = " + str(w) + ", n = " + str(n) + ", Iterations = "+str(Max_Iterations))
    plt.ylabel('y')
    plt.xlabel('x')
    plt.colorbar()
    plt.show()
