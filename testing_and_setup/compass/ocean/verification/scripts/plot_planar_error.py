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
nGrids = 3

ncfile = Dataset('output.nc', 'r')
#var = np.reshape(ncfile.variables['ssh'][0, :, iz], [ny, nx])

#dif = abs(var[1:ny - 1, 1:nx - 2] - sol[1:ny - 1, 1:nx - 2])
#err = abs((var[1:ny - 1, 1:nx - 2] - sol[1:ny - \
#          1, 1:nx - 2]) / sol[1:ny - 1, 1:nx - 2])
#difL2[i, j, k] = np.sqrt(np.mean(dif[:]**2))
#errL2[i, j, k] = np.sqrt(np.mean(err[:]**2))

fig = plt.gcf()
plt.clf()
fig.set_size_inches(20.0, 16.0)

#iTime=0
#xEdge = ncfile.variables['xEdge'][:]
#yEdge = ncfile.variables['yEdge'][:]
#angleEdge = ncfile.variables['angleEdge'][:]
#normalVelocity = ncfile.variables['normalVelocity'][0,:,0]
#print(np.size(normalVelocity))
#print(np.shape(normalVelocity))
#plt.scatter(xEdge,yEdge,s=80, c=normalVelocity, marker='.')
##plt.scatter(xEdge,yEdge,s=80, c=angleEdge, marker='.')
#plt.colorbar()
#plt.set_cmap('jet')

xCell = 1e-3*np.reshape(ncfile.variables['xCell'][:], [ny, nx])
yCell = 1e-3*np.reshape(ncfile.variables['yCell'][:], [ny, nx])

varNames = ['ssh','sshSolution']
varName='ssh'
ms=24
iTime=[0];
for j in range(len(iTime)):
      #for iVar, varName in enumerate(varNames):
      var = np.reshape(ncfile.variables[varName][iTime[j], :], [ny, nx])
      sol = np.reshape(ncfile.variables[varName+'Solution'][iTime[j], :], [ny, nx])
      dif = var-sol
      difL2 = np.sqrt(np.mean(dif[:]**2))
      difL1 = np.max(dif[:])
      print('ssh diffL2', difL2)
      print('ssh diffL1', difL1)
      #var_avg = var
      #sol_avg = sol
      #dif_avg = dif
      #for iy in range(0, ny, 2):
      #    for ix in range(1, nx - 2):
      #        var_avg[iy, ix] = (var[iy, ix + 1] + var[iy, ix]) / 2.0
      #        sol_avg[iy, ix] = (sol[iy, ix + 1] + sol[iy, ix]) / 2.0
      #        dif_avg[iy, ix] = (dif[iy, ix + 1] + dif[iy, ix]) / 2.0
      
      ax1 = plt.subplot(3, 3, 3 * j + 0 + 1)
      #plt.imshow(var_avg)
      plt.scatter(xCell,yCell,c=var,s=ms,marker='h')
      plt.title('model' + ' time = '+str(iTime[j]))
      plt.colorbar()
      plt.set_cmap('bwr')

      ax1 = plt.subplot(3, 3, 3 * j + 1 + 1)
      #plt.imshow(sol_avg)
      plt.scatter(xCell,yCell,c=sol,s=ms,marker='h')
      plt.title('sol' + ' time = '+str(iTime[j]))
      plt.colorbar()
      plt.set_cmap('bwr')

      ax1 = plt.subplot(3, 3, 3 * j + 2 + 1)
      #plt.imshow(dif_avg)
      plt.scatter(xCell,yCell,c=dif,s=ms,marker='h')
      plt.title('dif' + ' time = '+str(iTime[j]))
      plt.colorbar()
      plt.set_cmap('bwr')
#
#ax2 = plt.subplot(nGrids, 3, 3 * iTime + 2)
#plt.imshow(sol_avg[1:ny - 1, 1:nx - 2])
#plt.colorbar()
#plt.title(test + ' ' + tracer + ' ' + grid + ' solution')
#
#ax3 = plt.subplot(nGrids, 3, 3 * iTime + 3)
#plt.imshow(var_avg[1:ny - 1, 1:nx - 2] - sol_avg[1:ny - 1, 1:nx - 2])
#plt.colorbar()
#plt.title(test + ' ' + tracer + ' ' + grid + ' error')

ncfileMesh.close()
ncfile.close()

plt.savefig('plot.png')
