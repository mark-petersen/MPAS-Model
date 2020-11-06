#!/usr/bin/env python
'''
This script plots results from MPAS-Ocean convergence test.
'''
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


# mrp: read from file mesh after nx,ny attributes are added:
ncfileMesh = Dataset('planar_hex.nc', 'r')
nx = ncfileMesh.getncattr('nx')
ny = ncfileMesh.getncattr('ny')
iz = 1
iTime=0

ncfile = Dataset('output.nc', 'r')
var = np.reshape(ncfile.variables['ssh'][0, :], [ny, nx])
#var = np.reshape(ncfile.variables['ssh'][0, :, iz], [ny, nx])

#dif = abs(var[1:ny - 1, 1:nx - 2] - sol[1:ny - 1, 1:nx - 2])
#err = abs((var[1:ny - 1, 1:nx - 2] - sol[1:ny - \
#          1, 1:nx - 2]) / sol[1:ny - 1, 1:nx - 2])
#difL2[i, j, k] = np.sqrt(np.mean(dif[:]**2))
#errL2[i, j, k] = np.sqrt(np.mean(err[:]**2))

# --- Every other row in y needs to average two neighbors in x on planar hex mesh
var_avg = var
#sol_avg = sol
for iy in range(0, ny, 2):
    for ix in range(1, nx - 2):
        var_avg[iy, ix] = (var[iy, ix + 1] + var[iy, ix]) / 2.0
        #sol_avg[iy, ix] = (sol[iy, ix + 1] + sol[iy, ix]) / 2.0

#ax1 = plt.subplot(nGrids, 3, 3 * j + 1)
plt.imshow(var_avg)
plt.colorbar()
#plt.title(test + ' ' + tracer + ' ' + grid + ' computed')
#
#ax2 = plt.subplot(nGrids, 3, 3 * j + 2)
#plt.imshow(sol_avg[1:ny - 1, 1:nx - 2])
#plt.colorbar()
#plt.title(test + ' ' + tracer + ' ' + grid + ' solution')
#
#ax3 = plt.subplot(nGrids, 3, 3 * j + 3)
#plt.imshow(var_avg[1:ny - 1, 1:nx - 2] - sol_avg[1:ny - 1, 1:nx - 2])
#plt.colorbar()
#plt.title(test + ' ' + tracer + ' ' + grid + ' error')

ncfileMesh.close()
ncfile.close()

plt.savefig('plot.png')
