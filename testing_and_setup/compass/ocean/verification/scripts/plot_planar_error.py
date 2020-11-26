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
ncfileHex = Dataset('planar_hex.nc', 'r')
ncfileMesh = Dataset('initial_state.nc', 'r')
#ncfileBadAngleEdge = Dataset('mesh_bad_angleEdge.nc', 'r')
nx = ncfileHex.getncattr('nx')
ny = ncfileHex.getncattr('ny')
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
#normalVelocity = ncfile.variables['normalVelocity'][0,:,0]
#plt.scatter(xEdge,yEdge,s=80, c=normalVelocity, marker='.')
##plt.scatter(xEdge,yEdge,s=80, c=angleEdge, marker='.')
#plt.colorbar()
#plt.set_cmap('jet')

xCell = 1e-3*np.reshape(ncfileMesh.variables['xCell'][:], [ny, nx])
yCell = 1e-3*np.reshape(ncfileMesh.variables['yCell'][:], [ny, nx])

xEdge = 1e-3*ncfileMesh.variables['xEdge'][:]
yEdge = 1e-3*ncfileMesh.variables['yEdge'][:]
#angleEdge = 180./np.pi*ncfileBadAngleEdge.variables['angleEdge'][:]
angleEdge = 180./np.pi*ncfileMesh.variables['angleEdge'][:]

varNames = ['ssh','sshSolution']
varName='ssh'
varNameEdge='normalVelocity'
ms=24
ms=48
iTime=[0];
for i in range(len(iTime)):
      #for iVar, varName in enumerate(varNames):
      var = np.reshape(ncfile.variables[varName][iTime[i], :], [ny, nx])
      sol = np.reshape(ncfile.variables[varName+'Solution'][iTime[i], :], [ny, nx])
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
      
      j=0
      ax1 = plt.subplot(3, 3, 3 * j + 0 + 1)
      #plt.imshow(var_avg)
      plt.scatter(xCell,yCell,c=var,s=ms,marker='h')
      plt.title('model' + ' time = '+str(iTime[i]))
      plt.colorbar()
      plt.set_cmap('coolwarm')

      ax1 = plt.subplot(3, 3, 3 * j + 1 + 1)
      #plt.imshow(sol_avg)
      plt.scatter(xCell,yCell,c=sol,s=ms,marker='h')
      plt.title('sol' + ' time = '+str(iTime[i]))
      plt.colorbar()
      plt.set_cmap('coolwarm')

      ax1 = plt.subplot(3, 3, 3 * j + 2 + 1)
      #plt.imshow(dif_avg)
      plt.scatter(xCell,yCell,c=dif,s=ms,marker='h')
      plt.title('dif' + ' time = '+str(iTime[i]))
      plt.colorbar()
      plt.set_cmap('coolwarm')
# plot velocity on second row
      #for iVar, varName in enumerate(varNames):
      var = np.squeeze(ncfile.variables[varNameEdge][iTime[i], :, 1])
      sol = np.squeeze(ncfile.variables[varNameEdge+'Solution'][iTime[i], :, 1])
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
      
      j=1
      ax1 = plt.subplot(3, 3, 3 * j + 0 + 1)
      #plt.imshow(var_avg)
      plt.scatter(xEdge,yEdge,c=var,s=ms,marker='^')
      plt.title('model' + ' time = '+ varNameEdge + ' ' +str(iTime[i]))
      plt.colorbar()
      plt.set_cmap('coolwarm')

      ax1 = plt.subplot(3, 3, 3 * j + 1 + 1)
      #plt.imshow(sol_avg)
      plt.scatter(xEdge,yEdge,c=sol,s=ms,marker='^')
      plt.title('sol' + ' time = '+ varNameEdge + ' ' +str(iTime[i]))
      plt.colorbar()
      plt.set_cmap('coolwarm')

      ax1 = plt.subplot(3, 3, 3 * j + 2 + 1)
      #plt.imshow(dif_avg)
      plt.scatter(xEdge,yEdge,c=dif,s=ms,marker='^')
      plt.title('dif' + ' time = '+varNameEdge + ' ' +str(iTime[i]))
      plt.colorbar()
      plt.set_cmap('coolwarm')

      j=0
      plt.subplot(3, 3, 3 * j + 0 + 1)
      #plt.imshow(var_avg)
      plt.scatter(xEdge,yEdge,c=angleEdge,s=ms,marker='^')
      plt.title('model' + ' time = '+ varNameEdge + ' ' +str(iTime[i]))
      plt.colorbar()
      plt.set_cmap('hsv')
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
