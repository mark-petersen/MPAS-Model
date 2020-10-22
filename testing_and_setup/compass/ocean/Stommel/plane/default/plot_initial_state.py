#!/usr/bin/env python
'''
This script plots results from MPAS-Ocean convergence test.
'''
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

verbose=True
def comment(string):
    if verbose:
        print('***   '+string)

ncfile = Dataset('initial_state.nc', 'r')

comment('Convert all dimensions in netcdf file to local integer variables')
dimDict = ncfile.dimensions
for key in dimDict:
    vars()[key]= dimDict[key].size 
    if verbose:
        print('Adding dimension %s=%i'%(key,vars()[key]))


nCells=ncfile.dimensions['nCells'].size
nVertLevels=ncfile.dimensions['nVertLevels'].size
print('nCells',nCells)

comment('Find cartesian dimensions nx,ny')
try:
    history = ncfile.getncattr('history').split()
    nx = int(history[history.index('--nx')+1])
    ny = int(history[history.index('--ny')+1])
except:
    print("Error: Not able to read nx, ny from netcdf history attribute")
    print("history attribute: ",history)
    exit()
if nx*ny != nCells:
    print("Error: nx*ny != nCells. nx=%i, ny=%i, nx*ny=%i, nCells=%i"%(nx,ny,nx*ny,nCells))
    exit()
comment(' nCells = %i, nx = %i, ny = %i, nVertLevels = %i'%(nCells,nx,ny,nVertLevels))

#comment('Read in all variables')
ncVars = ncfile.variables
#for key in varDict:
#    var = varDict[key]
#    size = var.size
#    #print('key',key)
#    #print('var',var)
#    #print('var.size',var.size)
#    #print('var[:]',var[:])
#    #print('type(var',type(var))
#    #print('size',size)
#    print('var.dimensions',var.dimensions)
#    vars()[key]= varDict[key][:]
#    if verbose:
#        print('Adding variable %s'%(key))
#
#     var = np.reshape(ncfile.variables[tracer + 'Tend'][0, :, iz], [ny, nx])
#varList = ['xCell','yCell']
#xCell = np.reshape(ncfile.variables['xCell'][:], [ny, nx])
# --- Every other row in y needs to average two neighbors in x on planar hex mesh
#var_avg = var
#for iy in range(0, ny, 2):
#   for ix in range(1, nx - 2):
#       var_avg[iy, ix] = (var[iy, ix + 1] + var[iy, ix]) / 2.0
#print('xCell',xCell)
#print('np.size(xCell)',np.size(xCell))
#print('np.shape(xCell)',np.shape(xCell))
#            var = np.reshape(
#                ncfile.variables[tracer + 'Tend'][0, :, iz], [ny, nx])

fig = plt.gcf()
plt.clf()
fig.set_size_inches(20.0, 20.0)
iPlt=1



# 1D variables, function of x
varList = ['bottomDepth','maxLevelCell']
for varName in varList:
    plt.subplot(6, 3, iPlt)
    iPlt += 1
    plt.plot(ncVars['xCell'][:]*1e-3, ncVars[varName][:], '.')
    plt.title(varName)
    plt.xlabel('x location, km')
    plt.grid()
    plt.gca().invert_yaxis()

# 2D variables, function of (x,z)
varList = ['temperature','salinity','layerThickness']
for varName in varList:
    plt.subplot(6, 3, iPlt)
    iPlt += 1
    #plt.imshow(ncVars['xCell'][0:nx]*1e-3, ncVars[varName][0,0:nx,:])
    plt.imshow(ncVars[varName][0,0:nx,:].transpose())
    plt.title(varName)
    plt.xlabel('x location, km')
    plt.ylabel('Depth')
    plt.grid()
    plt.colorbar()
    plt.gca().invert_yaxis()
plt.savefig('variables.png')

#difL2 = np.zeros([nTests, nGrids, nTracers])
#errL2 = np.zeros([nTests, nGrids, nTracers])
#for i in range(nTests):
#    test = tests[i]
#
#    for k in range(nTracers):
#        fig = plt.gcf()
#        plt.clf()
#        fig.set_size_inches(20.0, 20.0)
#
#        for j in range(nGrids):
#            grid = grids[j]
#            ncfileIC = Dataset('../' + test + '_' + grid + '/init.nc', 'r')
#            ncfile = Dataset('../' + test + '_' + grid + '/output.nc', 'r')
#            tracer = tracers[k]
#            sol = np.reshape(
#                ncfileIC.variables[tracer + 'Tend'][0, :, iz], [ny, nx])
#            var = np.reshape(
#                ncfile.variables[tracer + 'Tend'][0, :, iz], [ny, nx])
#            dif = abs(var[1:ny - 1, 1:nx - 2] - sol[1:ny - 1, 1:nx - 2])
#            err = abs((var[1:ny - 1, 1:nx - 2] - sol[1:ny - \
#                      1, 1:nx - 2]) / sol[1:ny - 1, 1:nx - 2])
#            difL2[i, j, k] = np.sqrt(np.mean(dif[:]**2))
#            errL2[i, j, k] = np.sqrt(np.mean(err[:]**2))
#
#            # --- Every other row in y needs to average two neighbors in x on planar hex mesh
#            var_avg = var
#            sol_avg = sol
#            for iy in range(0, ny, 2):
#                for ix in range(1, nx - 2):
#                    var_avg[iy, ix] = (var[iy, ix + 1] + var[iy, ix]) / 2.0
#                    sol_avg[iy, ix] = (sol[iy, ix + 1] + sol[iy, ix]) / 2.0
#
#            ax1 = plt.subplot(nGrids, 3, 3 * j + 1)
#            plt.imshow(var_avg[1:ny - 1, 1:nx - 2])
#            plt.colorbar()
#            plt.title(test + ' ' + tracer + ' ' + grid + ' computed')
#            
#            ax2 = plt.subplot(nGrids, 3, 3 * j + 2)
#            plt.imshow(sol_avg[1:ny - 1, 1:nx - 2])
#            plt.colorbar()
#            plt.title(test + ' ' + tracer + ' ' + grid + ' solution')
#
#            ax3 = plt.subplot(nGrids, 3, 3 * j + 3)
#            plt.imshow(var_avg[1:ny - 1, 1:nx - 2] - sol_avg[1:ny - 1, 1:nx - 2])
#            plt.colorbar()
#            plt.title(test + ' ' + tracer + ' ' + grid + ' error')
#
#            ncfileIC.close()
#            ncfile.close()
#
#        plt.savefig(test + '_' + tracer + '_sections.png')
#
#
#fig = plt.gcf()
#plt.clf()
#fig.set_size_inches(20.0, 15.0)
#
#for i in range(nTests):
#    test = tests[i]
#    plt.subplot(nTests, 3, 3 * i + 1)
#    for k in range(len(tracers)):
#        tracer = tracers[k]
#        plt.loglog(dx, difL2[i, :, k], '-x', label=tracer)
#
#    plt.ylabel('diff: rms(exact[:] - calc[:])')
#    plt.legend()
#    plt.grid()
#plt.xlabel('cell width, km')
#
#for i in range(nTests):
#    test = tests[i]
#    plt.subplot(nTests, 3, 3 * i + 2)
#    for k in range(len(tracers)):
#        tracer = tracers[k]
#        plt.loglog(dx, errL2[i, :, k], '-x', label=tracer)
#
#    plt.title('Error in Redi tendancy term, ' + test)
#    plt.ylabel('error: rms((exact[:] - calc[:])/exact[:])')
#    plt.legend()
#    plt.grid()
#plt.xlabel('cell width, km')
#
#for i in range(nTests):
#    test = tests[i]
#    plt.subplot(nTests, 3, 3 * i + 3)
#    for k in range(len(tracers)):
#        tracer = tracers[k]
#        plt.loglog(dx, difL2[i, :, k] /
#                   np.max(difL2[i, :, k]), '-x', label=tracer)
#
#    plt.ylabel('normalized diff')
#    plt.grid()
#plt.xlabel('cell width, km')
#
#plt.savefig('convergence.png')
