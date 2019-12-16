#!/usr/bin/env python
'''
This script creates an initial condition file for MPAS-Ocean.
'''
from __future__ import absolute_import, division, print_function, \
    unicode_literals
import os
import shutil
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset


def main():

    shutil.copy2('base_mesh.nc','initial_state.nc')
    ds = Dataset('initial_state.nc', 'a', format='NETCDF3_64BIT_OFFSET')

    vertical_init(ds)
    tracer_init(ds)
    velocity_init(ds)
    coriolis_init(ds)
    others_init(ds)

    ds.close()


def vertical_init(ds):

    # config settings
    nVertLevels = 20
    thicknessAllLayers = 20 # [m] for evenly spaced layers

    # create new variables
    ds.createDimension('nVertLevels', nVertLevels)
    refLayerThickness = ds.createVariable('refLayerThickness', np.float64, ('nVertLevels',))
    maxLevelCell = ds.createVariable('maxLevelCell', np.int32, ('nCells',))
    refBottomDepth = ds.createVariable('refBottomDepth', np.float64, ('nVertLevels',))
    refZMid = ds.createVariable('refZMid', np.float64, ('nVertLevels',))
    bottomDepth = ds.createVariable('bottomDepth', np.float64, ('nCells',))
    bottomDepthObserved = ds.createVariable('bottomDepthObserved', np.float64, ('nCells',))
    layerThickness = ds.createVariable('layerThickness', np.float64, ('Time','nCells','nVertLevels',))
    restingThickness = ds.createVariable('restingThickness', np.float64, ('nCells','nVertLevels',))
    vertCoordMovementWeights = ds.createVariable('vertCoordMovementWeights', np.float64, ('nVertLevels',))

    # evenly spaced vertical grid
    refLayerThickness[:] = thicknessAllLayers

    # Create other variables from refLayerThickness
    refBottomDepth[0] = refLayerThickness[0]
    refZMid[0] = -0.5*refLayerThickness[0]
    for k in range(1,nVertLevels):
        refBottomDepth[k] = refBottomDepth[k-1] + refLayerThickness[k]
        refZMid[k] = -refBottomDepth[k-1]-0.5*refLayerThickness[k]
    vertCoordMovementWeights[:] = 1.0

    # flat bottom, no bathymetry
    maxLevelCell[:] = nVertLevels
    bottomDepth[:] = refBottomDepth[nVertLevels-1]
    bottomDepthObserved[:] = refBottomDepth[nVertLevels-1]
    for k in range(nVertLevels):
        layerThickness[0,:,k] = refLayerThickness[k]
        restingThickness[:,k] = refLayerThickness[k]


def tracer_init(ds):

    # config settings
    T0 = 20
    dTdx = 1e-4
    dTdy = 0
    dTdz = 10.0/1000

    S0 = 35
    dSdx = 0
    dSdy = 0
    dSdz = -10.0/1000

    # create new variables
    temperature = ds.createVariable('temperature', np.float64, ('Time','nCells','nVertLevels',))
    salinity = ds.createVariable('salinity', np.float64, ('Time','nCells','nVertLevels',))

    # obtain dimensions and mesh variables
    nVertLevels = len(ds.dimensions['nVertLevels'])
    nCells = len(ds.dimensions['nCells'])
    xCell = ds.variables['xCell']
    yCell = ds.variables['yCell']
    refZMid = ds.variables['refZMid']

    for iCell in range(0,nCells):
        for k in range(0,nVertLevels):
            temperature[0,iCell,k] \
                = T0 + xCell[iCell]*dTdx + yCell[iCell]*dTdy + refZMid[k]*dTdz
            salinity[0,iCell,k] \
                = S0 + xCell[iCell]*dSdx + yCell[iCell]*dSdy + refZMid[k]*dSdz


def velocity_init(ds):
    normalVelocity = ds.createVariable('normalVelocity', np.float64, ('Time','nEdges','nVertLevels',))
    normalVelocity = 0.0


def coriolis_init(ds):
    fEdge = ds.createVariable('fEdge', np.float64, ('nEdges',))
    fEdge[:] = 0.0
    fVertex = ds.createVariable('fVertex', np.float64, ('nVertices',))
    fVertex[:] = 0.0
    fCell = ds.createVariable('fCell', np.float64, ('nCells',))
    fCell[:] = 0.0


def others_init(ds):
    surfaceStress = ds.createVariable('surfaceStress', np.float64, ('Time','nEdges',))
    surfaceStress[:] = 0.0
    atmosphericPressure = ds.createVariable('atmosphericPressure', np.float64, ('Time','nCells',))
    atmosphericPressure[:] = 0.0
    boundaryLayerDepth = ds.createVariable('boundaryLayerDepth', np.float64, ('Time','nCells',))
    boundaryLayerDepth[:] = 0.0


if __name__ == '__main__':
    # If called as a primary module, run main
    main()
