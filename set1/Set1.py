import random
import numpy as np
import matplotlib.pyplot as plt


# Parameter settings
L = 1           # length of string
N = 100         # number of string parts
Dx = L/float(N) # Interval length
c = 1           # The constant C
Dt = 0.001      # Time step
time_limit = 1

# This is the matrix containing all the values of the string. So u[0] = u_i,0
u = np.zeros((int(time_limit/Dt)+1, N+1))           # sin(2pix)
v = np.zeros((int(time_limit/Dt)+1, N+1))           # sin(5pix)
w = np.zeros((int(time_limit/Dt)+1, N+1))           # The other one

# Initialisation
for i in range(1,N):
    u[0,i] = np.sin(2*np.pi*(Dx*i))
    v[0,i] = np.sin(5*np.pi*(Dx*i))
    if Dx*i > 0.2 and Dx*i < 0.4:
        w[0,i] = np.sin(5*np.pi*(Dx*i))
    else:
        w[0,i] = 0

# Keep track of the number of iterations and time.
time_step = 1
Time = 0
# Calculate the positions of the strings at different time values.
while Time < time_limit:
    for i in range(1,N):
        u[time_step, i] = (c*Dt/Dx)**2*(u[time_step-1, i+1] + u[time_step-1, i-1] - 2*u[time_step-1,i])-u[time_step-2, i]+2*u[time_step-1,i]
        v[time_step, i] = (c*Dt/Dx)**2*(v[time_step-1, i+1] + v[time_step-1, i-1] - 2*v[time_step-1,i])-v[time_step-2, i]+2*v[time_step-1,i]
        w[time_step, i] = (c*Dt/Dx)**2*(w[time_step-1, i+1] + w[time_step-1, i-1] - 2*w[time_step-1,i])-w[time_step-2, i]+2*w[time_step-1,i]

    Time += Dt
    time_step += 1


# Plotting sin(2pix)
plt.plot(u[3], label = "T = 0.003")
plt.plot(u[10], label = "T = 0.01")
plt.plot(u[50], label = "T = 0.05")
plt.plot(u[200], label = "T = 0.2")
plt.plot(u[1000], label = "T = 1")
plt.legend()
plt.title("Vibration amplitude at different times for $f(x) = sin(2 \pi x)$")
plt.xlabel("Position x")
plt.ylabel("Vibration Amplitude")
plt.show()

# plotting sin(5pix)
plt.plot(v[3], label = "T = 0.003")
plt.plot(v[10], label = "T = 0.01")
plt.plot(v[50], label = "T = 0.05")
plt.plot(v[200], label = "T = 0.2")
plt.plot(v[1000], label = "T = 1")
plt.legend()
plt.title("Vibration amplitude at different times for $g(x) sin(5 \pi x)$ ")
plt.xlabel("Position x")
plt.ylabel("Vibration Amplitude")
plt.show()

# plotting the other function
plt.plot(w[3], label = "T = 0.003")
plt.plot(w[10], label = "T = 0.01")
plt.plot(w[50], label = "T = 0.05")
plt.plot(w[200], label = "T = 0.2")
plt.plot(w[1000], label = "T = 1")
plt.legend()
plt.title("Vibration amplitude at different times for h(x)")
plt.xlabel("Position x")
plt.ylabel("Vibration Amplitude")
plt.show()
