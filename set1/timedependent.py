import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.special import erfc
import math

class Grid:
    def __init__(self, N, D, dt):
        self.N = N
        self.dx = 1.0/N
        self.dt = dt
        self.D = D
        self.const = self.dt*self.D/(self.dx**2)
        self.grid = np.zeros((N, N))

    def step(self):
        tempgrid = np.zeros((self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                value = self.get_boundry_values(i, j)
                tempgrid[i,j] = value

        self.grid = tempgrid

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

def analytical(x, D, t, i_limit=50):
    total = 0
    for i in range(i_limit):
        total += erfc((1-x+2*i)/(2*math.sqrt(D*t))) - erfc((1+x+2*i)/(2*math.sqrt(D*t)))

    return total

if __name__ == '__main__':
    N = 40
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

    analytic0001 = []
    analytic001 = []
    analytic01 = []
    analytic1 = []
    for i in range(N):
        analytic0001.append(analytical(i/float(N-1), D, 0.001))
        analytic001.append(analytical(i/float(N-1), D, 0.01))
        analytic01.append(analytical(i/float(N-1), D, 0.1))
        analytic1.append(analytical(i/float(N-1), D, 1.0))

    # analytic0001.append(1.0)
    # analytic001.append(1.0)
    # analytic01.append(1.0)
    # analytic1.append(1.0)
    x = [0, 5, 10, 15, 20, 25, 30, 35, 40]
    xticks = [i/float(N) for i in x]
    # plt.figure();
    # plt.plot(res[0], label="t = 0.001")
    # plt.plot(res[1], label="t = 0.01")
    # plt.plot(res[2], label="t = 0.1")
    # plt.plot(res[3], label="t = 1")
    # plt.legend()
    # plt.xlabel("y")
    # plt.ylabel("c")
    # plt.xticks(x, xticks, rotation='vertical')
    # plt.title("Numerical solution for concentration. t =" + `dt` + ", N = " + `N`)
    # plt.show()

    plt.figure()
    plt.plot(analytic0001, label="t=0.001")
    plt.plot(analytic001, label="t=0.01")
    plt.plot(analytic01, label="t=0.1")
    plt.plot(analytic1, label="t=1")
    plt.legend()
    plt.xlabel("y")
    plt.ylabel("c")
    plt.xticks(x, xticks, rotation='vertical')
    plt.title("Analytical solution for concentration. t =" + `dt` + ", N = " + `N`)
    plt.show()

    plt.figure()
    plt.plot(abs(analytic0001 - res[0]), label="absolute error t=0.001")
    plt.plot(abs(analytic001 - res[1]), label="absolute error t=0.01")
    plt.plot(abs(analytic01 - res[2]), label="absolute error t=0.1")
    plt.plot(abs(analytic1 - res[3]), label="absolute error t=1")
    plt.legend()
    plt.xlabel("y")
    plt.ylabel("c")
    plt.xticks(x, xticks, rotation='vertical')
    plt.title("Error analysis. t =" + `dt` + ", N = " + `N`)
    plt.show()
