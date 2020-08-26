import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

fig = plt.gcf()
iTime = [0]

# ---nx,ny, nz 
nx = 50
ny = 200
nz = 60


fig, ax = plt.subplots()

ncfile = Dataset('analysis_members/mpaso.hist.am.timeSeriesStatsMonthly.0001-01-01.nc', 'r')
var = ncfile.variables['timeMonthly_avg_velocityX']
var1 = np.reshape(var[1, :, :], [ny, nx, nz])
        # --- flip in y-dir
var = np.flipud(var1)
print('var=', var.shape)
print(var)

        # --- Every other row in y needs to average two neighbors in x on planar hex mesh
var_avg = var
for j in range(0, ny, 2):
   for i in range(1, nx-2):
     for k in range(0, nz-1):
        var_avg[j, i, k] = (var[j, i + 1, k] + var[j, i, k]) / 2.0

transect = np.reshape(var_avg[:, 20, :], [ny, nz]) 
dis = ax.imshow(transect.T)
ax.set_title("timeMonthly_avg_velocityX")
ax.set_xticks(np.arange(0, ny, step=25))
ax.set_yticks(np.arange(0, nz, step=10))

ax.set_xlabel('y, cells')
ax.set_ylabel('z, layers')
fig.colorbar(dis, ax=ax)
ncfile.close()

res = 20

plt.savefig("plots/velocityX_transect_timeSeriesStatsMonthly.0001-01-01.png")
