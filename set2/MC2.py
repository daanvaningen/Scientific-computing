import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
import math
# from numba import jit
from tqdm import tqdm

class Grid:
    def __init__(self, N, Epsilon, Max_Iterations, w=1, Object= 0, Pstick=1.0):
        self.MaxIt = Max_Iterations         # stop condition
        self.N = N                          # Grid size
        self.Iterations = 0                 # No. of current iteration
        self.Epsilon = Epsilon              # Epsilon (when do we accept change to be neglected.)
        self.grid = np.zeros((N, N))        # Making grid of n x n
        self.w = w
        self.const = w/4               # The w value for SOR iteration
        self.Object = Object                # A matrix containing the information of objects
        self.random_walkers = []
        self.init_random_walkers()
        self.candidates = []
        self.total_cluster = []
        self.get_object_edges(self.N/2, 0)
        self.stable = False
        self.Pstick = Pstick
        self.Hit = False

    # The grid will update untill certain criteria are met.
    def run(self):
        interval = 0
        for _ in tqdm(range(self.MaxIt)):
            self.move_random_walkers()

            if interval > 10:
                interval = 0
                tempgrid = np.zeros((self.N, self.N))
                for i in range(self.N):
                    for j in range(self.N):
                        value = self.get_boundry_values(i, j, tempgrid)
                        tempgrid[i,j] = value

                self.grid = tempgrid
            interval += 1

        return self.grid

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
        x = random.randint(0, self.N-1)
        y = 0
        return x,y

    def init_random_walkers(self):
        for _ in range(self.N):
            x, y = self.init_random_walker()
            self.random_walkers.append((x,y))

    def move_random_walkers(self):
        # Dont move up?
        for z in range(len(self.random_walkers)):
            rand_w = self.random_walkers[z]
            x_pos, y_pos = rand_w[0], rand_w[1]
            tries = 0
            while True:
                if(tries > 10):
                    x_pos, y_pos = self.init_random_walker()
                    self.random_walkers[z] = (x_pos, y_pos)
                    break
                x, y = self.rand_direction()
                if(not (x_pos + x, y_pos + y) in self.total_cluster):
                    x_pos += x
                    y_pos += y
                    break
                tries += 1


            # Check boundry
            if(y_pos >= self.N - 1 or y_pos < 0):
                x_pos, y_pos = self.init_random_walker()
                self.random_walkers[z] = (x_pos, y_pos)
            if(x_pos < 0):
                x_pos = self.N - 1
            elif(x_pos > self.N - 1):
                x_pos = 0

            self.random_walkers[z] = (x_pos, y_pos)
            # print self.rand_walk_x, self.rand_walk_y
            if (x_pos, y_pos) in self.candidates:
                # print self.rand_walk_x, self.rand_walk_y
                prob = random.random()
                if(prob < self.Pstick):
                    self.Object[x_pos, y_pos] = 1
                    nx, ny = self.init_random_walker()
                    self.random_walkers[z] = (nx, ny)
                    self.get_object_edges(x_pos, y_pos)


    def rand_direction(self):
        r = random.randint(0, 3)
        if r == 0:
            return 1, 0
        elif r == 1:
            return -1, 0
        elif r == 2:
            return 0,1
        else:
            return 0, -1

    def get_object_edges(self, x, y):
        self.total_cluster.append((x,y))
        neighbours = [(x+1, y), (x-1, y), (x, y+1), (x, y-1)]
        cluster = []
        for neighbour in neighbours:
            if(neighbour[0] < 0):
                neighbour = (self.N - 1, neighbour[1])
            elif(neighbour[0] > self.N - 1):
                neighbour = (0, neighbour[1])
            if(neighbour[1] < 0 or neighbour[1] > self.N - 1):
                continue

            if(self.Object[neighbour[0], neighbour[1]] == 1):
                cluster.append(neighbour)
            else:
                if not neighbour in self.candidates:
                    self.candidates.append(neighbour)

        for item in cluster:
            if not item in self.total_cluster:
                self.get_object_edges(item[0], item[1])


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
    N = 256
    w = 1.2
    Max_Iterations = 2000
    Grid_ = np.zeros((N,N))
    Grid_[N/2, 0] = 1
    saved = np.load('256x256_2.npy')
    MC = Grid(N, 10**-4, Max_Iterations, w, Object = Grid_, Pstick=1.0)
    MC.grid = saved
    MC.stable = True
    # print MC.candidates
    eq1 = MC.run()
    # print MC.total_cluster
    # eq1[0][MC.rand_walk_x, MC.rand_walk_y] = 1
    masked = np.ma.masked_where(MC.Object < 0.9, MC.Object)
    # #
    # x = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    # xticks = [i/float(N) for i in x]
    # yticks = [i/50 for i in range(50)]
    plt.imshow(np.rot90(eq1,1), origin = 'lower')
    plt.colorbar()
    plt.imshow(np.rot90(masked,1), cmap='Greys', interpolation=None)
    plt.title("MC DLA "+str(Max_Iterations)+' iterations')
    # plt.xticks(x, xticks)
    # plt.yticks(x, xticks)
    plt.ylabel('y')
    plt.xlabel('x')
    plt.show()

    # np.save('256x256_set2', Equilibrium3[0])
