#!/usr/bin/env python
'''
This script plots results from MPAS-Ocean convergence test.
'''
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import os
from scipy import interpolate
matplotlib.use('Agg')

varNames = ['ssh'] #,'normalVelocity']
#varNames = ['ssh','normalVelocity']
nVars = len(varNames)
plotDir = os.getcwd()
nRes = 5
nxGrid = 100

difL1 = np.zeros([nRes,nVars])
difL2 = np.zeros([nRes,nVars])

os.chdir('../../res1/simulation')
ncfileMesh = Dataset('initial_state.nc', 'r')
ncfileMesh.close()

dx = np.zeros([nRes])

iTime = 0
for iRes in range(nRes):
   os.chdir('../../res'+str(iRes+1)+'/simulation')
   ncfile = Dataset('output.nc', 'r')
   ncfileMesh = Dataset('initial_state.nc', 'r')
   xCell = ncfileMesh.variables['xCell'][:]
   yCell = ncfileMesh.variables['yCell'][:]
   #if iRes==0:
   #    x = np.arange(np.min(xCell),np.max(xCell),nGrid)
   #    y = np.arange(np.min(yCell),np.max(yCell),nGrid)
   #    xGrid, yGrid = np.meshgrid(x, y, sparse=True)

   dcEdge = 1e-3*ncfileMesh.variables['dcEdge'][:]
   dx[iRes] = np.max(dcEdge[:])

   for iVar, varName in enumerate(varNames):
       var = ncfile.variables[varName][:]
       if iRes==0:
           sol = ncfile.variables[varName+'Solution'][:]
           dif = var-sol
           xCell0 = xCell
           yCell0 = yCell
       else:
           fvar = interpolate.interp2d(xCell, yCell, var, kind='linear')
           dif = fvar(xCell0,yCell0) - sol

       difL1[iRes,iVar] = np.max(dif[:])
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
