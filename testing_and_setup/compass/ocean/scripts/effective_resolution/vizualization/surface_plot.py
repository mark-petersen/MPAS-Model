import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

fig = plt.gcf()
iTime = [0]

# ---nx,ny 
nx = 50
ny = 200
ncfile = Dataset('analysis_members/mpaso.hist.am.timeSeriesStatsMonthly.0001-01-01.nc', 'r')
var = ncfile.variables['timeMonthly_avg_velocityX']


fig, ax = plt.subplots()

var1 = np.reshape(var[1, :, 0], [ny, nx])
        # --- flip in y-dir
var = np.flipud(var1)

        # --- Every other row in y needs to average two neighbors in x on planar hex mesh
var_avg = var
for j in range(0, ny, 2):
    for i in range(1, nx - 2):
        var_avg[j, i] = (var[j, i + 1] + var[j, i]) / 2.0

dis = ax.imshow(
    var_avg)
ax.set_title("timeMonthly_avg_velocityX")
ax.set_xticks(np.arange(0, nx, step=25))
ax.set_yticks(np.arange(0, ny, step=25))

ax.set_xlabel('x, cells')
ax.set_ylabel('y, cells')
fig.colorbar(dis, ax=ax)
ncfile.close()

res = 20

plt.savefig("plots/velocityX_surface_timeSeriesStatsMonthly.0001-01-01.png")
