import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

L = 1           # length of string
N = 100         # number of string parts
Dx = L/float(N) # Interval length
c = 1           # The constant C
Dt = 0.01      # Time step
time_limit = 1

# This is the matrix containing all the values of the string. So u[0] = u_i,0
cur = np.zeros(N)           # sin(2pix)
nex = np.zeros(N)

for i in range(1,N):
    cur[i] = np.sin(2*np.pi*(Dx*i))

fig = plt.figure()
ax = plt.axes(xlim=(0, 100), ylim=(-2, 2))
# ax = plt.axes()
line, = ax.plot([], [], lw=2)

def init():
    global cur
    x = np.linspace(0, N, N)
    line.set_data(x, cur)
    return line,

# animation function.  This is called sequentially
def animate(i):
    global cur
    global nex
    x = np.linspace(0, N, N)

    for j in range(1,N):
        if(j == 0 or j == N - 1):
            nex[j] = 0
        else:
            nex[j] = (c*Dt/Dx)**2*(cur[j+1] + cur[j-1] - 2*cur[j])-cur[j] + 2 * cur[j]

    line.set_data(x, nex)
    cur = nex
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=20, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.xlabel('x')
plt.ylabel('y')
plt.show()
