import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
import math
from numba import jit

class Grid:
    def __init__(self, N, Epsilon, Max_Iterations, w=1, Object= 0):
        self.MaxIt = Max_Iterations         # stop condition
        self.N = N                          # Grid size
        self.Iterations = 0                 # No. of current iteration
        self.Epsilon = Epsilon              # Epsilon (when do we accept change to be neglected.)
        self.grid = np.zeros((N, N))        # Making grid of n x n
        self.w = w
        self.const = w/4               # The w value for SOR iteration
        self.deltas = []                    # The convergence measures
        self.Object = Object                # A matrix containing the information of objects
        self.init_random_walker()

    # The grid will update untill certain criteria are met.
    def run(self):
        Stable = False
        while not Stable and self.Iterations < self.MaxIt:
            print self.Iterations
            self.move_random_walker()
            self.Iterations += 1
            tempgrid = np.zeros((self.N, self.N))
            Stable = True           # Unless determined otherwise
            delta = 0               # check the largest convergence measure
            for i in range(self.N):
                for j in range(self.N):
                    value = self.get_boundry_values(i, j, tempgrid)
                    tempgrid[i,j] = value
                    if abs(value - self.grid[i,j]) > self.Epsilon:  # Means not yet stable
                        Stable = False
                        if abs(value - self.grid[i,j]) > delta:
                            delta = abs(value - self.grid[i,j])

            self.deltas.append(delta)
            self.grid = tempgrid
        return self.grid, self.deltas

    # @jit
    def get_boundry_values(self, i, j, tempgrid):
        # Checking whether it is part of an object
        if self.Object[i,j] == 1:
            return 0
        else:
            # The values when we do a jacobi iteration
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

    def init_random_walker(self):
        self.rand_walk_x = random.randint(0, self.N-1)
        self.rand_walk_y = 0

    def move_random_walker(self):
        r = random.randint(0.3):
        if r == 0:
            self.rand_walk_x += 1
        elif r == 1:
            self.rand_walk_x -= 1
        elif r == 2:
            self.rand_walk_y += 1
        elif r == 3:
            self.rand_walk_y -= 1

        # Check boundry
        if(self.rand_walk_y >= self.N - 1 or self.rand_walk_y < 0):
            self.init_random_walker()
        elif(self.rand_walk_x < 0):
            self.rand_walk_x = self.N - 1
        elif(self.rand_walk_x > self.N - 1):
            self.rand_walk_x = 0

        
def ObjectCreation(size):
    grid = np.zeros((N, N))
    radius = int(size/2.0)
    middle = int(N/2.0)

    for i in range(middle-radius, middle+radius):
        for j in range(middle-radius, middle+radius):
            grid[i,j] = 1

    return grid

if __name__ == '__main__':
    # Initial conditions
    N = 50
    w = 1.2
    Max_Iterations = 100
    Grid_ = np.zeros((N,N))
    Grid_[N/2, N-1] = 1
    MC = Grid(N, 0, Max_Iterations, w, Object = Grid_)
    eq1 = MC.run()
    eq1[0][MC.rand_walk_x, MC.rand_walk_y] = 1

    x = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    xticks = [i/float(N) for i in x]
    yticks = [i/50 for i in range(50)]
    plt.imshow(np.rot90(eq1[0],1), origin = 'lower')
    plt.title("MC DLA")
    # plt.xticks(x, xticks)
    # plt.yticks(x, xticks)
    plt.ylabel('y')
    plt.xlabel('x')
    plt.colorbar()
    plt.show()
