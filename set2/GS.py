import numpy as np
import random
import matplotlib.pyplot as plt
import math
from tqdm import tqdm

class Grid:
    def __init__(self, N, dt, dx, Du, Dv, f, k):
        self.N = N
        self.dt = dt
        self.dx = dx
        self.Du = Du
        self.Dv = Dv
        self.f = f
        self.k = k
        self.constU = self.dt*self.Du/self.dx**2
        self.constV = self.dt*self.Dv/self.dx**2
        self.grid = np.empty((self.N,self.N, 2))

    def init_grid(self, size, noise=0.001, rand_size=5):
        size = size + random.randint(-2,2)
        radius = int(size/2.0)
        middle = int(self.N/2.0)
        for i in range(N):
            for j in range(N):
                if(i > middle + radius or i < middle - radius):
                    self.grid[i,j,0], self.grid[i,j,1] = 0.5, 0.0
                else:
                    if(j > middle + radius or j < middle - radius):
                        self.grid[i,j,0], self.grid[i,j,1] = 0.5, 0.0
                    else:
                        self.grid[i,j,0], self.grid[i,j,1] = 0.5, 0.25

    def step(self):
        tempgrid = np.empty((self.N,self.N, 2))
        for i in range(self.N):
            for j in range(self.N):
                tempgrid[i,j,0] = self.get_boundry_values(i, j, 'u')
                tempgrid[i,j,1] = self.get_boundry_values(i, j, 'v')

        self.grid = tempgrid

    def get_boundry_values(self, i, j, c_type):
        if(c_type == 'u'):
            if(i == 0 or j == 1 or i == self.N - 1 or j == self.N - 1):
                return 0.5
            else:
                return self.grid[i, j, 0] + self.constU*(self.grid[i+1, j, 0]+
                    self.grid[i-1, j, 0] + self.grid[i, j+1, 0] + self.grid[i, j-1, 0] \
                    - 4*self.grid[i,j,0]) - self.grid[i,j,0]*self.grid[i,j,1]**2 \
                    + self.f*(1 - self.grid[i,j,0])

        elif(c_type == 'v'):
            if(i == 0 or j == 1 or i == self.N - 1 or j == self.N - 1):
                return 0.0
            else:
                return self.grid[i, j, 1] + self.constV*(self.grid[i+1, j, 1]+
                    self.grid[i-1, j, 1] + self.grid[i, j+1, 1] + self.grid[i, j-1, 1] \
                    - 4*self.grid[i,j,1]) + self.grid[i,j,0]*self.grid[i,j,1]**2 \
                    - self.grid[i,j,1]*(self.f + self.k)

if __name__ == '__main__':
    N = 100
    dt = 1.0
    dx = 1.0
    Du = 0.16 # 0.16
    Dv = 0.08 # 0.08
    f = 0.035 # 0.035
    k = 0.1 # 0.060
    GS = Grid(N, dt, dx, Du, Dv, f, k)
    GS.init_grid(10, 1)

    iterations = 5000
    for _ in tqdm(range(iterations)):
        GS.step()

    u = np.zeros((N, N))
    v = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            u[i,j] = GS.grid[i,j][0]
            v[i,j] = GS.grid[i,j][1]

    plt.imshow(u)
    plt.imshow(v)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Gray Scott diffusion, '+str(iterations)+' iterations \n $D_u = '+str(Du)+'$,\
 $D_v = '+str(Dv)+'$, $f = '+str(f)+'$, $k = '+str(k)+'$' )
    plt.show()
