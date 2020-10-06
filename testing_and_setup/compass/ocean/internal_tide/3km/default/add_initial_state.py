#!/usr/bin/env python
'''
This script creates an initial condition file for MPAS-Ocean.

Internal tide test case
see:
Demange, J., Debreu, L., Marchesiello, P., Lemarié, F., Blayo, E., Eldred, C., 2019. Stability analysis of split-explicit free surface ocean models: Implication of the depth-independent barotropic mode approximation. Journal of Computational Physics 398, 108875. https://doi.org/10.1016/j.jcp.2019.108875
Marsaleix, P., Auclair, F., Floor, J.W., Herrmann, M.J., Estournel, C., Pairaud, I., Ulses, C., 2008. Energy conservation issues in sigma-coordinate free-surface ocean models. Ocean Modelling 20, 61–89. https://doi.org/10.1016/j.ocemod.2007.07.005
'''
import os
import shutil
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import argparse
import math
import time
verbose=True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', dest='input_file',
                        default='base_mesh.nc',
                        help='Input file, containing base mesh'
                        )
    parser.add_argument('-o', '--output_file', dest='output_file',
                        default='initial_state.nc',
                        help='Output file, containing initial variables'
                        )
    parser.add_argument('-L', '--nVertLevels', dest='nVertLevels',
                        default=50,
                        help='Number of vertical levels'
                        )
    nVertLevels = parser.parse_args().nVertLevels

    input_file = parser.parse_args().input_file
    output_file = parser.parse_args().output_file
    shutil.copy2(input_file, output_file)
    ds = Dataset(output_file, 'a', format='NETCDF3_64BIT_OFFSET')

    maxDepth = 5000.0 

    comment('obtain dimensions and mesh variables')
    nCells = ds.dimensions['nCells'].size 
    xCell = ds.variables['xCell']
    xEdge = ds.variables['xEdge']
    xVertex = ds.variables['xVertex']
    yCell = ds.variables['yCell']
    yEdge = ds.variables['yEdge']
    yVertex = ds.variables['yVertex']

# for some reason this didn't work:
    # obtain dimensions and mesh variables
    #for key in ds.dimensions:
    #    vars()[key]= ds.dimensions[key].size 
    #for key in ds.variables:
    #    vars()[key]= ds.variables[key]
    #print('xEdge',xEdge[:])


    comment('create new variables')
    ds.createDimension('nVertLevels', nVertLevels)
    refLayerThickness = ds.createVariable(
        'refLayerThickness', np.float64, ('nVertLevels',))
    maxLevelCell = ds.createVariable('maxLevelCell', np.int32, ('nCells',))
    refBottomDepth = ds.createVariable(
        'refBottomDepth', np.float64, ('nVertLevels',))
    refZMid = ds.createVariable('refZMid', np.float64, ('nVertLevels',))
    bottomDepth = ds.createVariable('bottomDepth', np.float64, ('nCells',))
    ssh = ds.createVariable('ssh', np.float64, ('nCells',))
    bottomDepthObserved = ds.createVariable(
        'bottomDepthObserved', np.float64, ('nCells',))
    layerThickness = ds.createVariable('layerThickness', np.float64, ('Time', 'nCells', 'nVertLevels',))
    zMid = ds.createVariable('zMid', np.float64, ('Time', 'nCells', 'nVertLevels',))
    restingThickness = ds.createVariable(
        'restingThickness', np.float64, ('nCells', 'nVertLevels',))
    vertCoordMovementWeights = ds.createVariable(
        'vertCoordMovementWeights', np.float64, ('nVertLevels',))

    # Adjust coordinates so first edge is at zero in x and y
    xOffset = min(xEdge)
    xCell -= xOffset
    xEdge -= xOffset
    xVertex -= xOffset
    yOffset = min(yEdge)
    yCell -= yOffset
    yEdge -= yOffset
    yVertex -= yOffset

    # For periodic domains, the max cell coordinate is also the domain width
    Lx = max(xCell)
    xMid = min(xCell) + 0.5*(max(xCell) - min(xCell))

    comment('create reference variables for z-level grid')
    refLayerThickness[:] = maxDepth/nVertLevels
    refBottomDepth[0] = refLayerThickness[0]
    refZMid[0] = -0.5 * refLayerThickness[0]
    for k in range(1, nVertLevels):
        refBottomDepth[k] = refBottomDepth[k - 1] + refLayerThickness[k]
        refZMid[k] = -refBottomDepth[k - 1] - 0.5 * refLayerThickness[k]
    comment('z-level: ssh in top layer only')
    vertCoordMovementWeights[:] = 0.0
    vertCoordMovementWeights[0] = 1.0

    for iCell in range(0, nCells):
        x = xCell[iCell]
        # Marsaleix et al 2008 page 81
        # Gaussian function in depth for deep sea ridge
        bottomDepth[iCell] = 5000.0 - 1000.0*math.exp( -(( x - xMid)/150e3)**2 )
        # SSH varies from 0 to 1m across the domain
        ssh[iCell] = xCell[iCell]/4800e3
        #print('iCell',iCell)

        # z-star: spread layer thicknesses proportionally
        #layerThickness[0, iCell, :] = refLayerThickness[:]*(maxDepth+ssh[iCell])/maxDepth

        # z-level: ssh in top layer only
        layerThickness[0, iCell, :] = refLayerThickness[:]
        layerThickness[0, iCell, 0] += ssh[iCell]

        for k in range(nVertLevels-1,0,-1):
            # kF is the Fortran index, starts at 1.
            kF = k+1
            if bottomDepth[iCell] > refBottomDepth[k-1]:
                maxLevelCell[iCell] = kF
                # Partial bottom cells
                layerThickness[0, iCell, k] = bottomDepth[iCell] - refBottomDepth[k-1]
                zMid[0, iCell, k] = -bottomDepth[iCell] + 0.5*layerThickness[0, iCell, k]
                break
            else:
                layerThickness[0, iCell, k] = -1.0
                zMid[0, iCell, k] = 0.0

        for k in range(maxLevelCell[iCell]-2,-1,-1):
            zMid[0, iCell, k] = zMid[0, iCell, k+1]  \
               + 0.5*(layerThickness[0, iCell, k+1] + layerThickness[0, iCell, k])

    restingThickness[:, :] = layerThickness[0, :, :]
    bottomDepthObserved[:] = bottomDepth[:]


    comment('initialize tracers')
    rho0 = 1000.0 # kg/m^3
    rhoz = -2.0e-4 # kg/m^3/m in z 
    S0 = 35.0

#test for now:
    time1 = time.time()
    salinity = ds.createVariable('salinity', np.float64, ('Time', 'nCells', 'nVertLevels',))
    print(' time 0: %f'%((time.time()-time1)))

    time1 = time.time()
    for iCell in range(0, nCells):
        for k in range(0, maxLevelCell[iCell]):
            salinity[0, iCell, k] = S0
    print(' time 1: %f'%((time.time()-time1)))

    time1 = time.time()
    salinity[:] = S0
    print(' time 2: %f'%((time.time()-time1)))
#test for now end


    # linear equation of state
    # rho = rho0 - alpha*(T-Tref) + beta*(S-Sref)
    # set S=Sref
    # T = Tref - (rho - rhoRef)/alpha 
    config_eos_linear_alpha = 0.2
    config_eos_linear_beta = 0.8
    config_eos_linear_Tref = 10.0
    config_eos_linear_Sref = 35.0
    config_eos_linear_densityref = 1000.0

    # create new variables
    temperature = ds.createVariable('temperature', np.float64, ('Time', 'nCells', 'nVertLevels',))
    salinity = ds.createVariable('salinity', np.float64, ('Time', 'nCells', 'nVertLevels',))
    density = ds.createVariable('density', np.float64, ('Time', 'nCells', 'nVertLevels',))

    for iCell in range(0, nCells):
        for k in range(0, maxLevelCell[iCell]):
            salinity[0, iCell, k] = S0
            density[0,iCell,k] = rho0 + rhoz*zMid[0, iCell, k]
            # T = Tref - (rho - rhoRef)/alpha 
            temperature[0,iCell,k] = config_eos_linear_Tref \
                - (density[0,iCell,k] - config_eos_linear_densityref)/config_eos_linear_alpha

    comment('initialize velocity')
    normalVelocity = ds.createVariable('normalVelocity', np.float64, ('Time', 'nEdges', 'nVertLevels',))
    normalVelocity[:] = 0.0

    comment('initialize coriolis terms')
    fEdge = ds.createVariable('fEdge', np.float64, ('nEdges',))
    fEdge[:] = 0.0
    fVertex = ds.createVariable('fVertex', np.float64, ('nVertices',))
    fVertex[:] = 0.0
    fCell = ds.createVariable('fCell', np.float64, ('nCells',))
    fCell[:] = 0.0

    comment('initialize other fields')
    surfaceStress = ds.createVariable('surfaceStress', np.float64, ('Time', 'nEdges',))
    surfaceStress[:] = 0.0
    atmosphericPressure = ds.createVariable('atmosphericPressure', np.float64, ('Time', 'nCells',))
    atmosphericPressure[:] = 0.0
    boundaryLayerDepth = ds.createVariable('boundaryLayerDepth', np.float64, ('Time', 'nCells',))
    boundaryLayerDepth[:] = 0.0

    comment('close file')
    ds.close()


def comment(string):
    if verbose:
        print('***   '+string)

if __name__ == '__main__':
    # If called as a primary module, run main
    main()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
