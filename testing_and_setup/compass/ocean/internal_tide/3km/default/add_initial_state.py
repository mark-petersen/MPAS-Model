#!/usr/bin/env python
'''
This script creates an initial condition file for MPAS-Ocean.

Internal tide test case
see:
Demange, J., Debreu, L., Marchesiello, P., Lemarié, F., Blayo, E., Eldred, C., 2019. Stability analysis of split-explicit free surface ocean models: Implication of the depth-independent barotropic mode approximation. Journal of Computational Physics 398, 108875. https://doi.org/10.1016/j.jcp.2019.108875
Marsaleix, P., Auclair, F., Floor, J.W., Herrmann, M.J., Estournel, C., Pairaud, I., Ulses, C., 2008. Energy conservation issues in sigma-coordinate free-surface ocean models. Ocean Modelling 20, 61–89. https://doi.org/10.1016/j.ocemod.2007.07.005
'''
# import packages {{{
import os
import shutil
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import argparse
import math
# }}}


def main(): # {{{ 
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

    vertical_init(ds, nVertLevels)
    tracer_init(ds)
    velocity_init(ds)
    coriolis_init(ds)
    others_init(ds)

    ds.close()
# }}}


def vertical_init(ds, nVertLevels):
    maxDepth = 5000.0 
    # {{{

    # create new variables
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
    layerThickness = ds.createVariable(
        'layerThickness', np.float64, ('Time', 'nCells', 'nVertLevels',))
    restingThickness = ds.createVariable(
        'restingThickness', np.float64, ('nCells', 'nVertLevels',))
    vertCoordMovementWeights = ds.createVariable(
        'vertCoordMovementWeights', np.float64, ('nVertLevels',))

    # obtain dimensions and mesh variables
    nCells = len(ds.dimensions['nCells'])
    xCell = ds.variables['xCell']
    xEdge = ds.variables['xEdge']
    xVertex = ds.variables['xVertex']
    yCell = ds.variables['yCell']
    yEdge = ds.variables['yEdge']
    yVertex = ds.variables['yVertex']

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

    # create reference variables for z-level grid
    refLayerThickness[:] = maxDepth/nVertLevels
    refBottomDepth[0] = refLayerThickness[0]
    refZMid[0] = -0.5 * refLayerThickness[0]
    for k in range(1, nVertLevels):
        refBottomDepth[k] = refBottomDepth[k - 1] + refLayerThickness[k]
        refZMid[k] = -refBottomDepth[k - 1] - 0.5 * refLayerThickness[k]
    vertCoordMovementWeights[:] = 0.0
    vertCoordMovementWeights[0:int(nVertLevels/2)] = 1.0

    for iCell in range(0, nCells):
        x = xCell[iCell]
        # Marsaleix et al 2008 page 81
        # Gaussian function in depth for deep sea ridge
        bottomDepth[iCell] = 5000.0 - 1000.0*math.exp( -(( x - xMid)/15e3)**2 )
        # SSH varies from 0 to 1m across the domain
        ssh[iCell] = xCell[iCell]/4800e3
        for k in range(nVertLevels,0,-1):
            if bottomDepth[iCell] > refBottomDepth[k-1]:
                # Convert python zero-based to Fortran by adding 1
                maxLevelCell[iCell] = k+1
                layerThickness[0, iCell, k] = bottomDepth[iCell] - refBottomDepth[k-1]
                break
            else:
                layerThickness[0, iCell, k] = -99
        # spread layer thicknesses proportionally, with z-star
        layerThickness[0, iCell, :] = refLayerThickness[:]*(maxDepth+ssh[iCell])/maxDepth

    restingThickness[:, :] = layerThickness[0, :, :]
    bottomDepthObserved[:] = bottomDepth[:]

# }}}


def tracer_init(ds, thicknessAllLayers):
    rho0 = 1000.0 # kg/m^3
    rhoz = -2.0e-4 # kg/m^3/m in z 
    S0 = 35.0
# {{{
    h = thicknessAllLayers
    # create new variables
    temperature = ds.createVariable(
        'temperature', np.float64, ('Time', 'nCells', 'nVertLevels',))
    salinity = ds.createVariable(
        'salinity', np.float64, ('Time', 'nCells', 'nVertLevels',))
    layerThickness = ds.variables['layerThickness']

    # obtain dimensions and mesh variables # {{{
    nVertLevels = len(ds.dimensions['nVertLevels'])
    nCells = len(ds.dimensions['nCells'])
    xCell = ds.variables['xCell']
    yCell = ds.variables['yCell']
    # For periodic domains, the max cell coordinate is also the domain width
    Lx = max(xCell)
    Ly = max(yCell)
    refZMid = ds.variables['refZMid']
    refBottomDepth = ds.variables['refBottomDepth']
    H = max(refBottomDepth)

    for iCell in range(0, nCells):
        x = xCell[iCell]
        y = yCell[iCell]
        for k in range(0, nVertLevels):
            z = refZMid[k]

            salinity[0, iCell, k] = S0
            temperature[0,iCell,k] = Tx*(x + x0) + Ty + Tz*z
# }}}

def velocity_init(ds):
    # {{{
    normalVelocity = ds.createVariable(
        'normalVelocity', np.float64, ('Time', 'nEdges', 'nVertLevels',))
    normalVelocity[:] = 0.0
# }}}


def coriolis_init(ds):
    # {{{
    fEdge = ds.createVariable('fEdge', np.float64, ('nEdges',))
    fEdge[:] = 0.0
    fVertex = ds.createVariable('fVertex', np.float64, ('nVertices',))
    fVertex[:] = 0.0
    fCell = ds.createVariable('fCell', np.float64, ('nCells',))
    fCell[:] = 0.0
# }}}


def others_init(ds):
    # {{{
    surfaceStress = ds.createVariable(
        'surfaceStress', np.float64, ('Time', 'nEdges',))
    surfaceStress[:] = 0.0
    atmosphericPressure = ds.createVariable(
        'atmosphericPressure', np.float64, ('Time', 'nCells',))
    atmosphericPressure[:] = 0.0
    boundaryLayerDepth = ds.createVariable(
        'boundaryLayerDepth', np.float64, ('Time', 'nCells',))
    boundaryLayerDepth[:] = 0.0
# }}}


if __name__ == '__main__':
    # If called as a primary module, run main
    main()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
