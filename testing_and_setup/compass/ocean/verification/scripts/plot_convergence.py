#!/usr/bin/env python
'''
This script plots results from MPAS-Ocean convergence test.
'''
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import os
from scipy.interpolate import RegularGridInterpolator

matplotlib.use('Agg')

varNames = ['ssh'] #,'normalVelocity']
#varNames = ['ssh','normalVelocity']
nVars = len(varNames)
plotDir = os.getcwd()
nRes = 5
nGrid = 10

difL1 = np.zeros([nRes,nVars])
difL2 = np.zeros([nRes,nVars])

os.chdir('../../res1/simulation')
#ncfileMesh = Dataset('initial_state.nc', 'r')
#ncfileMesh.close()

dx = np.zeros([nRes])

iTime = 0
for iRes in range(nRes):
   os.chdir('../../res'+str(iRes+1)+'/simulation')
   ncfileHex = Dataset('planar_hex.nc', 'r')
   nx = ncfileHex.getncattr('nx')
   ny = ncfileHex.getncattr('ny')
   ncfileHex.close()
   ncfile = Dataset('output.nc', 'r')
   ncfileMesh = Dataset('initial_state.nc', 'r')
   xCell = np.reshape( ncfileMesh.variables['xCell'][:], [ny, nx])
   yCell = np.reshape( ncfileMesh.variables['yCell'][:], [ny, nx])
   yCell1D = yCell[:,0]
   xCell1D = np.zeros(2*nx-2)
   for ix in range(0, nx-1):
       xCell1D[2*ix] = (xCell[0, ix] + xCell[0, ix+1]) / 2.0
       xCell1D[2*ix+1] = xCell[0, ix+1]

   if iRes==0:
       YGridRef, XGridRef = np.meshgrid(yCell1D[1:nx-2], xCell1D[1:nx-2])
       #x = np.arange(np.min(xCell),np.max(xCell),nGrid)
       #y = np.arange(np.min(yCell),np.max(yCell),nGrid)
       #xGrid, yGrid = np.meshgrid(x, y, sparse=True)

   dcEdge = 1e-3*ncfileMesh.variables['dcEdge'][:]
   dx[iRes] = np.max(dcEdge[:])

   for iVar, varName in enumerate(varNames):
       var = np.reshape( ncfile.variables[varName][iTime, :], [ny, nx])
       sol = np.reshape( ncfile.variables[varName+'Solution'][iTime, :], [ny, nx])
       # --- Every other row in y needs to average two neighbors in x on planar hex mesh
       varGrid = np.zeros([ny,2*nx-2])
       solGrid = np.zeros([ny,2*nx-2])
       for iy in range(0, ny, 2):
           for ix in range(0, nx-1):
               varGrid[iy, 2*ix] = (var[iy, ix] + var[iy, ix+1]) / 2.0
               varGrid[iy, 2*ix+1] = var[iy, ix+1]
               varGrid[iy+1, 2*ix] = var[iy+1, ix]
               varGrid[iy+1, 2*ix+1] = (var[iy+1, ix] + var[iy+1, ix+1]) / 2.0
               solGrid[iy, 2*ix] = (sol[iy, ix] + sol[iy, ix+1]) / 2.0
               solGrid[iy, 2*ix+1] = sol[iy, ix+1]
               solGrid[iy+1, 2*ix] = sol[iy+1, ix]
               solGrid[iy+1, 2*ix+1] = (sol[iy+1, ix] + sol[iy+1, ix+1]) / 2.0
       if iRes==0:
           dif = varGrid-solGrid
       else:
           varFxn = RegularGridInterpolator((yCell1D, xCell1D), varGrid)
           solFxn = RegularGridInterpolator((yCell1D, xCell1D), solGrid)
           varInterp = varFxn((YGridRef, XGridRef))
           solInterp = solFxn((YGridRef, XGridRef))
           #solInterp = solFxn((yCell1Dref, xCell1Dref))
           dif = varInterp - solInterp

       difL1[iRes,iVar] = np.max(abs(dif[:]))
       difL2[iRes,iVar] = np.sqrt(np.mean(dif[:]**2))
         
   ncfile.close()
   ncfileMesh.close()

os.chdir(plotDir)
print('ssh diffL1', difL1)
print('ssh diffL2', difL2)

print('dx[i+1]/dx[i], (dx[i+1]/dx[i])^2, difL1[i+1]/difL1[i], difL2[i+1]/difL2[i]')
for i in range(nRes-1):
    print(dx[i+1]/dx[i], (dx[i+1]/dx[i])**2, difL1[i+1]/difL1[i], difL2[i+1]/difL2[i])

# ax1 = plt.subplot(3, 3, 3 * j + 0 + 1)
# #plt.imshow(var_avg)
# plt.scatter(xCell,yCell,c=var,s=ms,marker='h')
# plt.title('model' + ' time = '+str(iTime[j]))
# plt.colorbar()
# plt.set_cmap('bwr')
# 
# ax1 = plt.subplot(3, 3, 3 * j + 1 + 1)
# #plt.imshow(sol_avg)
# plt.scatter(xCell,yCell,c=sol,s=ms,marker='h')
# plt.title('sol' + ' time = '+str(iTime[j]))
# plt.colorbar()
# plt.set_cmap('bwr')
# 
# ax1 = plt.subplot(3, 3, 3 * j + 2 + 1)
# #plt.imshow(dif_avg)
# plt.scatter(xCell,yCell,c=dif,s=ms,marker='h')
# plt.title('dif' + ' time = '+str(iTime[j]))
# plt.colorbar()
# plt.set_cmap('bwr')
# 
# plt.savefig('plot.png')
