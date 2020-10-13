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
import xarray as xr
from mpas_tools.io import write_netcdf
import argparse
import math
import time
verbose = True


def main():
    timeStart = time.time()
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
    ds = xr.open_dataset(parser.parse_args().input_file)

    #comment('obtain dimensions and mesh variables')
    nCells = ds['nCells'].size
    nEdges = ds['nEdges'].size
    nVertices = ds['nVertices'].size

    xCell = ds['xCell']
    xEdge = ds['xEdge']
    xVertex = ds['xVertex']
    yCell = ds['yCell']
    yEdge = ds['yEdge']
    yVertex = ds['yVertex']

    # Adjust coordinates so first edge is at zero in x and y
    xOffset = xEdge.min()
    xCell -= xOffset
    xEdge -= xOffset
    xVertex -= xOffset
    yOffset = np.min(yEdge)
    yCell -= yOffset
    yEdge -= yOffset
    yVertex -= yOffset

    comment('create and initialize variables')
    time1 = time.time()

    varsZ = [ 'refLayerThickness', 'refBottomDepth', 'refZMid', 'vertCoordMovementWeights']
    for var in varsZ:
        globals()[var] = np.nan * np.ones(nVertLevels)

    vars2D = ['ssh', 'bottomDepth', 'bottomDepthObserved',
        'surfaceStress', 'atmosphericPressure', 'boundaryLayerDepth']
    for var in vars2D:
        globals()[var] = np.nan * np.ones(nCells)
    maxLevelCell = np.ones(nCells, dtype=np.int32)

    vars3D = [
        'temperature', 'salinity',
        'layerThickness', 'restingThickness', 'zMid',
        'density']
    for var in vars3D:
        globals()[var] = np.nan * np.ones([1, nCells, nVertLevels])

    # reference vertical grid spacing
    maxDepth = 5000.0
    refLayerThickness[:] = maxDepth / nVertLevels
    refBottomDepth[0] = refLayerThickness[0]
    refZMid[0] = -0.5 * refLayerThickness[0]
    for k in range(1, nVertLevels):
        refBottomDepth[k] = refBottomDepth[k - 1] + refLayerThickness[k]
        refZMid[k] = -refBottomDepth[k - 1] - 0.5 * refLayerThickness[k]

    # Marsaleix et al 2008 page 81
    # Gaussian function in depth for deep sea ridge
    xMid = 0.5 * (min(xCell) + max(xCell))
    bottomDepth[:] = 5000.0 - 1000.0 * np.exp(-((xCell[:] - xMid) / 150e3)**2)
    # SSH varies from 0 to 1m across the domain
    ssh[:] = xCell[:] / 4800e3

    # Compute maxLevelCell and layerThickness for z-level (variation only on top)
    vertCoordMovementWeights[:] = 0.0
    vertCoordMovementWeights[0] = 1.0
    for iCell in range(0, nCells):
        for k in range(nVertLevels - 1, 0, -1):
            if bottomDepth[iCell] > refBottomDepth[k - 1]:
                maxLevelCell[iCell] = k
                # Partial bottom cells
                layerThickness[0, iCell, k] = bottomDepth[iCell] - refBottomDepth[k - 1]
                break
        layerThickness[0, iCell, 0:maxLevelCell[iCell] ] = refLayerThickness[0:maxLevelCell[iCell]]
        layerThickness[0, iCell, 0] += ssh[iCell]

    # Compute zMid (same, regardless of vertical coordinate)
    for iCell in range(0, nCells):
        k = maxLevelCell[iCell]
        zMid[0, iCell, k] = -bottomDepth[iCell] + \
            0.5 * layerThickness[0, iCell, k]
        for k in range(maxLevelCell[iCell] - 1, -1, -1):
            zMid[0, iCell, k] = zMid[0, iCell, k + 1] + 0.5 * \
                (layerThickness[0, iCell, k + 1] + layerThickness[0, iCell, k])
    restingThickness[:, :] = layerThickness[0, :, :]
    restingThickness[:, 0] = refLayerThickness[0]
    bottomDepthObserved[:] = bottomDepth[:]

    # initialize tracers
    rho0 = 1000.0  # kg/m^3
    rhoz = -2.0e-4  # kg/m^3/m in z
    S0 = 35.0

    # linear equation of state
    # rho = rho0 - alpha*(T-Tref) + beta*(S-Sref)
    # set S=Sref
    # T = Tref - (rho - rhoRef)/alpha
    config_eos_linear_alpha = 0.2
    config_eos_linear_beta = 0.8
    config_eos_linear_Tref = 10.0
    config_eos_linear_Sref = 35.0
    config_eos_linear_densityref = 1000.0

    for k in range(0, nVertLevels):
        activeCells = k <= maxLevelCell
        salinity[0, activeCells, k] = S0
        density[0, activeCells, k] = rho0 + rhoz * zMid[0, activeCells, k]
        # T = Tref - (rho - rhoRef)/alpha
        temperature[0, activeCells, k] = config_eos_linear_Tref \
            - (density[0, activeCells, k] - config_eos_linear_densityref) / \
              config_eos_linear_alpha

    # initial velocity on edges
    normalVelocity = (('Time', 'nEdges', 'nVertLevels',), 0.0)

    # Coriolis parameter
    ds['fCell'] = (('nCells', 'nVertLevels',), np.zeros([nCells, nVertLevels]))
    ds['fEdge'] = (('nEdges', 'nVertLevels',), np.zeros([nEdges, nVertLevels]))
    ds['fVertex'] = (('nVertices', 'nVertLevels',), np.zeros([nVertices, nVertLevels]))

    # surface fields
    surfaceStress[:] = 0.0
    atmosphericPressure[:] = 0.0
    boundaryLayerDepth[:] = 0.0
    print('   time: %f' % ((time.time() - time1)))

    comment('finalize and write file')
    time1 = time.time()
    ds['maxLevelCell'] = (['nCells'], maxLevelCell + 1)
    for var in varsZ:
        ds[var] = (['nVertLevels'], globals()[var])
    for var in vars2D:
        ds[var] = (['nCells'], globals()[var])
    for var in vars3D:
        ds[var] = (['Time', 'nCells', 'nVertLevels'], globals()[var])
    # If you prefer not to have NaN as the fill value, you should consider
    # using mpas_tools.io.write_netcdf() instead
    ds.to_netcdf('initial_state.nc', format='NETCDF3_64BIT_OFFSET')
    # write_netcdf(ds,'initial_state.nc')
    print('   time: %f' % ((time.time() - time1)))
    print('Total time: %f' % ((time.time() - timeStart)))


def comment(string):
    if verbose:
        print('***   ' + string)


if __name__ == '__main__':
    # If called as a primary module, run main
    main()
