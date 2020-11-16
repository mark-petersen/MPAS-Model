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
test = 'coastal Kelvin Wave, '
time_step_method = 'split_explicit'
time_step_method = 'RK4'
#varNames = ['ssh','normalVelocity']
nVars = len(varNames)
plotDir = os.getcwd()
nRes = 5

difL1 = np.zeros([nRes,nVars])
difL2 = np.zeros([nRes,nVars])
dx = np.zeros([nRes])

iTime = 0
for iRes in range(nRes):
   os.chdir('../../resolution'+str(iRes+1)+'/simulation_'+ time_step_method)
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
# Choose the grid for comparison.
# Compare only where data exists for the exact solution
       xCell1DRef = xCell[0,1:nx-2]
       yCell1DRef = yCell[2:ny-2:2,0]
       YGridRef, XGridRef = np.meshgrid(yCell1DRef, xCell1DRef)

   dcEdge = 1e-3*ncfileMesh.variables['dcEdge'][:]
   dx[iRes] = np.max(dcEdge[:])

   for iVar, varName in enumerate(varNames):
       var = np.reshape( ncfile.variables[varName][iTime, :], [ny, nx])
       sol = np.reshape( ncfile.variables[varName+'Solution'][iTime, :], [ny, nx])
       error_on_same_grid = True

       if (error_on_same_grid):
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
               dif = varInterp - solInterp
       else:
           dif = var - sol

       difL1[iRes,iVar] = np.max(abs(dif[:]))
       difL2[iRes,iVar] = np.sqrt(np.mean(dif[:]**2))
         
   ncfile.close()
   ncfileMesh.close()

os.chdir(plotDir)
print('ssh diffL1', difL1)
print('ssh diffL2', difL2)

print('dx[i+1]/dx[i], (dx[i+1]/dx[i])^2, difL1[i+1]/difL1[i], difL2[i+1]/difL2[i], L1Slope, L2Slope')
for i in range(nRes-1):
    print(dx[i+1]/dx[i], (dx[i+1]/dx[i])**2, difL1[i+1]/difL1[i], difL2[i+1]/difL2[i], difL1[i+1]/difL1[i]/(dx[i+1]/dx[i]), difL2[i+1]/difL2[i]/(dx[i+1]/dx[i]))

p = np.polyfit(np.log10(dx),np.log10(difL2),1)
conv = abs(p[0])

yfit = dx**p[0]*10**p[1]

fig, ax = plt.subplots()
plt.loglog(dx,yfit,'k')
plt.loglog(dx,difL2,'or')
plt.annotate('Order of Convergence = {}'.format(np.round(conv,2)),xycoords='axes fraction',xy=(0.04,0.95),fontsize=14)
plt.xlabel('cell width [km]',fontsize=14)
plt.ylabel('L2 Norm',fontsize=14)
plt.title(test + time_step_method)
#ax.grid()
plt.savefig('convergence_'+test+time_step_method+'.png',bbox_inches='tight', pad_inches=0.1)
