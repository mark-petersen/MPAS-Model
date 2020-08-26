import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from matplotlib import animation

fig = plt.gcf()

# ---nx,ny 
nx = 50
ny = 200
#ncfile = Dataset('output.nc', 'r')
#var = ncfile.variables['kineticEnergyCell']
#print('variable_size', var.shape)

fig, ax = plt.subplots()
#ax.set_title("kineticEnergyCell")
ax.set_xticks(np.arange(0, nx, step=10))
ax.set_yticks(np.arange(0, ny, step=10))
ax.set_xlabel('x, cells')
ax.set_ylabel('y, cells')
#line, = ax.plot([], [], lw=2)

#def init():
#    line.set_data([], [])
#    return line,

#var_test = []
#for k in range(0,29):
def animate(k):
    ncfile = Dataset('output.nc', 'r')
    var = ncfile.variables['kineticEnergyCell']

    var1 = np.reshape(var[k, :, 0], [ny, nx])
        # --- flip in y-dir
    var = np.flipud(var1)

        # --- Every other row in y needs to average two neighbors in x on planar hex mesh
    var_avg = var
    for j in range(0, ny, 2):
        for i in range(1, nx - 2):
            var_avg[j, i] = (var[j, i + 1] + var[j, i]) / 2.0

    dis = ax.imshow(var_avg)
    ax.set_title("kineticEnergyCell"+ str(k))
    ncfile.close()
    #lineS.set_data(var_avg)
    #return lineS,
    return dis,

anim = animation.FuncAnimation(fig, animate, frames=30)

#dis = ax.imshow(lineS)
#fig.colorbar(dis, ax=ax)
#ncfile.close()

#anim = animation.FuncAnimation(fig, animate, frames=30)



#plt.show()
anim.save("plots/animation_kineticEnergyCell.mp4")
