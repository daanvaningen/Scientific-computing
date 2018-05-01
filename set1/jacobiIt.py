import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
import math

class Grid:
    def __init__(self, N, Epsilon, Max_Iterations, form, w=1, Object= 0):
        self.MaxIt = Max_Iterations         # stop condition
        self.N = N                          # Grid size
        self.Iterations = 0                 # No. of current iteration
        self.Epsilon = Epsilon              # Epsilon (when do we accept change to be neglected.)
        self.grid = np.zeros((N, N))        # Making grid of n x n
        self.w = w                          # The w value for SOR iteration
        self.deltas = []                    # The convergence measures
        self.Object = Object                # A matrix containing the information of objects
        
        if form == "SOR":
            self.const = w/4.0              # Value dependent on W
        else:
            self.const = 1/4.0              # Constant value for jacobi Iterations

        self.form  = form

    # The grid will update untill certain criteria are met.
    def run(self):
        Stable = False
        while not Stable and self.Iterations < self.MaxIt:
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

    def get_boundry_values(self, i, j, tempgrid):
        # Checking whether it is part of an object
        if self.Object[i,j] ==1:
            return 0
        else:
            # The values when we do a jacobi iteration
            if self.form == "Gauss" or self.form == "SOR":
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

            # The values when we do jacobi iterations
            elif self.form == "Jac":
                if(j == 0 or j == self.N-1):
                    return j/(self.N-1)
                
                elif(i == 0):
                    return self.grid[0, j] + self.const*(self.grid[1, j] +
                            self.grid[N-1, j] + self.grid[0, j+1] + self.grid[0, j-1] - 4*self.grid[0,j])
                elif(i == self.N-1):
                    return self.grid[i, j] + self.const*(self.grid[1, j] +
                            self.grid[i-1, j] + self.grid[i, j+1] + self.grid[i, j-1] - 4*self.grid[i,j])
                else:
                    return self.grid[i, j] + self.const*(self.grid[i+1, j] +
                            self.grid[i-1, j] + self.grid[i, j+1] + self.grid[i, j-1] - 4*self.grid[i,j])


def ObjectCreation(size):
    grid = np.zeros((N, N))
    radius = int(size/2.0)
    middle = int(N/2.0)

    for i in range(middle-radius, middle+radius):
        for j in range(middle-radius, middle+radius):
            grid[i,j] = 1

    return grid

# The visuals of the jacobi iteration.
def J():
    # Grid creation for jacobi
    Jacobi10_3 = Grid(N, 10**(-3), Max_Iterations, "Jac")
    Jacobi10_5 = Grid(N, 10**(-3.5), Max_Iterations, "Jac")
    Jacobi10_7 = Grid(N, 10**(-4), Max_Iterations, "Jac")

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
    plt.title("Jacobi Iteration")
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
    plt.colorbar()
    plt.show()


# The visuals of the Gauss iteration.
def G():
    # Grid creation for jacobi
    Gauss10_3 = Grid(N, 10**(-3), Max_Iterations, "Gauss")
    Gauss10_5 = Grid(N, 10**(-3.5), Max_Iterations, "Gauss")
    Gauss10_7 = Grid(N, 10**(-4), Max_Iterations, "Gauss")

    # Run the model
    Equilibrium3 = Gauss10_3.run()
    Equilibrium5 = Gauss10_5.run()
    Equilibrium7 = Gauss10_7.run()

    # Iterations it took to obtain epsilon.
    epsilon10_3 = Equilibrium3[0][0]   # One vertical line representing changes in y
    epsilon10_5 = Equilibrium5[0][0]
    epsilon10_7 = Equilibrium7[0][0]

    x = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    xticks = [i/float(N) for i in x]

    # Show the change in y
    yticks = [i/50 for i in range(50)]
    plt.title("Gauss Iteration")
    plt.plot(yticks, epsilon10_3, label = "Epsilon = 10^-3")
    plt.plot(yticks, epsilon10_5, label = "Epsilon = 10^-35")
    plt.plot(yticks, epsilon10_7, label = "Epsilon = 10^-4")
    plt.legend()
    plt.xlabel("y")
    plt.ylabel("Concentration c")
    plt.show()

    plt.imshow(np.rot90(Equilibrium3[0],3), origin = 'lower')
    plt.title("Stable state Gauss iteration with epsilon = 10^-3")
    plt.xticks(x, xticks)
    plt.yticks(x, xticks)
    plt.ylabel('y')
    plt.xlabel('x')
    plt.colorbar()
    plt.show()

def S():
    # Grid creation for jacobi
    SOR10_3 = Grid(N, 0, Max_Iterations, "SOR", w-0.5, Object = Grid_)
    SOR10_5 = Grid(N, 0, Max_Iterations, "SOR", w, Object = Grid_)
    SOR10_7 = Grid(N, 0, Max_Iterations, "SOR", w+0.5, Object = Grid_)

    # Run the model
    Equilibrium3 = SOR10_3.run()
    Equilibrium5 = SOR10_5.run()
    Equilibrium7 = SOR10_7.run()

    epsilon10_3 = Equilibrium3[0][0]   # One vertical line representing changes in y
    epsilon10_5 = Equilibrium5[0][0]
    epsilon10_7 = Equilibrium7[0][0]

    x = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    xticks = [i/float(N) for i in x]

    # Show the change in y
    yticks = [i/50 for i in range(50)]
    plt.title("SOR iteration with max "+str(Max_Iterations)+ " iterations")
    plt.plot(yticks, epsilon10_3, label = "w =" +str(round(w-0.5,1)))
    plt.plot(yticks, epsilon10_5, label = "w =" +str(w))
    plt.plot(yticks, epsilon10_7, label = "w =" +str(w+0.5))
    plt.legend()
    plt.xlabel("y")
    plt.ylabel("Concentration c")
    plt.show()

    plt.imshow(np.rot90(Equilibrium7[0],3), origin = 'lower')
    plt.title("SOR with object")
    plt.xticks(x, xticks)
    plt.yticks(x, xticks)
    plt.ylabel('y')
    plt.xlabel('x')
    plt.colorbar()
    plt.show()
    
def Compare():
    SOR_ = Grid(N, 0, Max_Iterations, "SOR", w)
    GAUSS_ = Grid(N, 0, Max_Iterations, "Gauss")
    JACOBI_ = Grid(N, 0, Max_Iterations, "Jac")

    # Run the model
    Equilibrium3 = SOR_.run()
    Equilibrium5 = GAUSS_.run()
    Equilibrium7 = JACOBI_.run()

    # One vertical line representing changes in y
    Sor_change = Equilibrium3[0][0]
    Gauss_change = Equilibrium5[0][0]
    Jacobi_change = Equilibrium7[0][0]

    x = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    xticks = [i/float(N) for i in x]

    # Show the change in y
    yticks = [i/50 for i in range(50)]
    plt.title("Comparance with max " +str(Max_Iterations)+" iterations")
    plt.plot(yticks, Sor_change, label = "SOR with w= " + str(w))
    plt.plot(yticks, Gauss_change, label = "Gauss")
    plt.plot(yticks, Jacobi_change, label = "Jacobi")
    plt.legend()
    plt.xlabel("y")
    plt.ylabel("Concentration c")
    plt.show()

def Convergence_Measure(function = "Jac", w=1):
    if function == "Jac":
        jac = Grid(N, 0, Max_Iterations, "Jac")
        data1 = jac.run()
        delta_values1 = data1[1]

        gauss = Grid(N, 0, Max_Iterations, "Gauss")
        data2 = gauss.run()
        delta_values2 = data2[1]
        
        iterations = [x for x in range(1, len(delta_values1)+1)]

        plt.plot(iterations, delta_values1, label = "Jacobi")
        plt.plot(iterations, delta_values2, label = "Gauss")
        plt.title("Jacobi and Gauss convergence")
        plt.legend()
        plt.xlabel("Number of iterations k")
        plt.ylabel("$\delta$")
        plt.yscale('log')
        plt.show()

    else:
        SOR1 = Grid(N, 0, Max_Iterations, "SOR", 0.5)
        data1 = SOR1.run()
        delta_values1 = data1[1]

        SOR2 = Grid(N, 0, Max_Iterations, "SOR", 1.0)
        data2 = SOR2.run()
        delta_values2 = data2[1]

        SOR3 = Grid(N, 0, Max_Iterations, "SOR", 1.5)
        data3 = SOR3.run()
        delta_values3 = data3[1]
        
        iterations = [x for x in range(1, len(delta_values1)+1)]

        plt.plot(iterations, delta_values1, label = "w = 0.5")
        plt.plot(iterations, delta_values2, label = "w = 1.0")
        plt.plot(iterations, delta_values3, label = "w = 1.5")
        plt.title("SOR convergence for different w")
        plt.legend()
        plt.xlabel("Number of iterations k")
        plt.ylabel("$\delta$")
        plt.yscale('log')
        plt.show()


def Error():
    w = 1.7
    Best_error = 10000
    Square_errors = []
    Ws = []
    best_w = 1.7
    while w <= 1.94:
        SOR = Grid(N, 0, 300, "SOR", w, Grid_)
        Data = SOR.run()
        c_values = Data[0][0]               # The estimated C values
        yticks = [i/N for i in range(N)]    # The analytical C values

        # error calculation
        error = 0
        for y in yticks:
            error += (y-c_values[int(y*N)])**2

        if error < Best_error:
            Best_error = error
            best_w = w

        Ws.append(w)
        Square_errors.append(error)
        
        w+=0.01

    print("Best w found was " +str(best_w) + "with an error of: ", Best_error)

    plt.plot(Ws, Square_errors)
    plt.xlabel("W")
    plt.ylabel("Square error")
    plt.yscale('log')
    plt.title("Error measured for different W")
    plt.show()
    
if __name__ == '__main__':
    # Initial conditions
    N = 50
    w = 1.2                       # For the SOR model
    Max_Iterations = 500
    Max_Iterations_SOR = 500        # Max iterations specifically for SOR
    Jacobi = False                  # True if you want results, false if u want to skip
    Gauss = False                   # True if you want results, false if u want to skip
    SOR = True                    # True if you want results, false if u want to skip
    Comparance = False
    Convergence = False             
    ErrorCalculation = True
    ObjectGrid = True

    if ObjectGrid:
        Grid_ = ObjectCreation(10)
    else:
        Grid_ = np.zeros((N, N))
        
    if Jacobi:
        J()

    if Gauss:
        G()

    if SOR:
        S()

    if Comparance:
        Compare()

    if Convergence:
        Convergence_Measure()
        Convergence_Measure(function = "SOR")

    if ErrorCalculation:
        Error()
