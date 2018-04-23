import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
import math

class Grid:
    def __init__(self, N, Epsilon, Max_Iterations):
        self.MaxIt = Max_Iterations         # stop condition
        self.N = N                          # Grid size
        self.Iterations = 0                 # No. of current iteration
        self.Epsilon = Epsilon              # Epsilon (when do we accept change to be neglected.)
        self.grid = np.zeros((N, N))        # Making grid of n x n
        self.const = 1/4.0                  # Constant value for jacobi Iterations

    # The grid will update untill certain criteria are met.
    def run(self):
        Stable = False
        while not Stable and self.Iterations < self.MaxIt:
            self.Iterations += 1
            tempgrid = np.zeros((self.N, self.N))
            Stable = True           # Unless determined otherwise
            for i in range(self.N):
                for j in range(self.N):
                    value = self.get_boundry_values(i, j)
                    tempgrid[i,j] = value
                    if abs(value - self.grid[i,j]) > self.Epsilon:  # Means not yet stable
                        Stable = False

            self.grid = tempgrid
        return self.grid, self.Iterations

    def get_boundry_values(self, i, j):
        if(j == 0 or j == self.N-1):
            return j/(self.N-1)
        elif(i == 0):
            return self.grid[0, j] + self.const*(self.grid[1, j] +
                    self.grid[N-1, j] + self.grid[0, j+1] + self.grid[0, j-1] - 4*self.grid[0,j])
        elif(i == self.N-1):
            return self.grid[i, j] + self.const*(self.grid[0, j] +
                    self.grid[i-1, j] + self.grid[i, j+1] + self.grid[i, j-1] - 4*self.grid[i,j])
        else:
            return self.grid[i, j] + self.const*(self.grid[i+1, j] +
                    self.grid[i-1, j] + self.grid[i, j+1] + self.grid[i, j-1] - 4*self.grid[i,j])

if __name__ == '__main__':
    # Initial conditions
    N = 50
    Max_Iterations = 1000000000

    # Grid creation
    Jacobi10_3 = Grid(N, 10**(-3), Max_Iterations)
    Jacobi10_5 = Grid(N, 10**(-3.5), Max_Iterations)
    Jacobi10_7 = Grid(N, 10**(-4), Max_Iterations)

    # Run the model
    Equilibrium3 = Jacobi10_3.run()
    Equilibrium5 = Jacobi10_5.run()
    Equilibrium7 = Jacobi10_7.run()

    # Iterations it took to obtain epsilon.
    epsilon10_3 = Equilibrium3[0][0]   # One vertical line representing changes in y
    epsilon10_5 = Equilibrium5[0][0]
    epsilon10_7 = Equilibrium7[0][0]

    x = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    xticks = [i/float(N) for i in x]

    # Show the change in y
    yticks = [i/50 for i in range(50)]
    plt.plot(yticks, epsilon10_3, label = "Epsilon = 10^-3")
    plt.plot(yticks, epsilon10_5, label = "Epsilon = 10^-35")
    plt.plot(yticks, epsilon10_7, label = "Epsilon = 10^-4")
    plt.legend()
    plt.xlabel("y")
    plt.ylabel("Concentration c")
    plt.show()

    plt.imshow(np.rot90(Equilibrium3[0],3), origin = 'lower')
    plt.title("Stable state with Jacobi iteration with epsilon = 10^-3")
    plt.xticks(x, xticks)
    plt.yticks(x, xticks)
    plt.ylabel('y')
    plt.xlabel('x')
    plt.legend()
    plt.show()
