import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.special import erfc

class Grid:
    def __init__(self, N, D, dt):
        self.N = N
        self.dx = 1.0/N
        self.dt = dt
        self.D = D
        self.const = self.dt*self.D/(self.dx**2)
        self.grid = np.zeros((N, N))
        print(self.const)

    def step(self):
        tempgrid = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                value = self.get_boundry_values(i, j)
                tempgrid[i,j] = value

        self.grid = tempgrid

    def get_boundry_values(self, i, j):
        if(j == 0 or j == N-1):
            return j/(N-1)
        elif(i == 0):
            return self.grid[0, j] + self.const*(self.grid[1, j] +
                    self.grid[N-1, j] + self.grid[0, j+1] + self.grid[0, j-1] - 4*self.grid[0,j])
        elif(i == N-1):
            return self.grid[i, j] + self.const*(self.grid[0, j] +
                    self.grid[i-1, j] + self.grid[i, j+1] + self.grid[i, j-1] - 4*self.grid[i,j])
        else:
            return self.grid[i, j] + self.const*(self.grid[i+1, j] +
                    self.grid[i-1, j] + self.grid[i, j+1] + self.grid[i, j-1] - 4*self.grid[i,j])

if __name__ == '__main__':
    N = 30
    D = 1
    dt = 0.0001
    grid = Grid(N, D, dt)
    T = 0
    t0001 = True
    t001 = True
    t01 = True
    t1 = True

    res = []

    while(T < 1.0):
        grid.step()
        T += dt
        if(t0001 and T > 0.001):
            t0001 = False
            res.append(grid.grid[0])
        elif(t001 and T > 0.01):
            t001 = False
            res.append(grid.grid[0])
        elif(t01 and T > 0.1):
            t01 = False
            res.append(grid.grid[0])

    res.append(grid.grid[0])

    plt.figure();
    plt.plot(res[0], label="t = 0.001")
    plt.plot(res[1], label="t = 0.01")
    plt.plot(res[2], label="t = 0.1")
    plt.plot(res[3], label="t = 1")
    plt.legend()
    plt.xlabel("y")
    plt.ylabel("c")
    plt.title("Numerical solution for concentration. t =" + `dt` + ", N = " + `N`)
    plt.show()
