#!/usr/bin/env python
'''
This script creates the analytical solution file (Stommel's test case)
'''
# import packages {{{
import os
import shutil
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import argparse
# }}}


def main():
    # {{{
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', dest='input_file',
                        default='base_mesh.nc',
                        help='Input file, containing base mesh'
                        )
    parser.add_argument('-o', '--output_file', dest='output_file',
                        default='velx.nc',
                        help='Output file, containing initial variables'
                        )
    parser.add_argument('-L', '--nVertLevels', dest='nVertLevels',
                        default=1,
                        type=int,
                        help='Number of vertical levels'
                        )
    parser.add_argument('-H',
        '--thicknessAllLayers',
        dest='thicknessAllLayers',
        default=200,
        type=float,
        help='thickness of each layer, [m]')
    thicknessAllLayers = parser.parse_args().thicknessAllLayers
    nVertLevels = parser.parse_args().nVertLevels

    input_file = parser.parse_args().input_file
    output_file = parser.parse_args().output_file
    shutil.copy2(input_file, output_file)
    ds = Dataset(output_file, 'a', format='NETCDF3_64BIT_OFFSET')

    vertical_init(ds, thicknessAllLayers, nVertLevels)

    ds.close()
# }}}


def vertical_init(ds, thicknessAllLayers, nVertLevels):
    # {{{

    # create new variables # {{{
    ds.createDimension('nVertLevels', nVertLevels)
    refLayerThickness = ds.createVariable(
        'refLayerThickness', np.float64, ('nVertLevels',))
    maxLevelCell = ds.createVariable('maxLevelCell', np.int32, ('nCells',))
    refBottomDepth = ds.createVariable(
        'refBottomDepth', np.float64, ('nVertLevels',))
    refZMid = ds.createVariable('refZMid', np.float64, ('nVertLevels',))
    bottomDepth = ds.createVariable('bottomDepth', np.float64, ('nCells',))
    bottomDepthObserved = ds.createVariable(
        'bottomDepthObserved', np.float64, ('nCells',))
    layerThickness = ds.createVariable(
        'layerThickness', np.float64, ('Time', 'nCells', 'nVertLevels',))
    restingThickness = ds.createVariable(
        'restingThickness', np.float64, ('nCells', 'nVertLevels',))
    vertCoordMovementWeights = ds.createVariable(
        'vertCoordMovementWeights', np.float64, ('nVertLevels',))
    windStressZonal = ds.createVariable(
        'windStressZonal', np.float64, ('Time', 'nCells',))
    velx = ds.createVariable(
            'velx', np.float64, ('nCells',))
    vely = ds.createVariable(
                'vely', np.float64, ('nCells',))
    sshx = ds.createVariable(
                'sshx', np.float64, ('nCells',))


    # }}}

    # evenly spaced vertical grid
    refLayerThickness[:] = thicknessAllLayers
    # make first layer deep to avoid z^2 derivative problems near zero.
    refLayerThickness[0] = 200
    
    nVertLevels = len(ds.dimensions['nVertLevels'])
    nCells = len(ds.dimensions['nCells'])
    lonCell = ds.variables['lonCell']
    latCell = ds.variables['latCell']
    # For periodic domains, the max cell coordinate is also the domain width
    Lx = max(lonCell)
    Ly = max(latCell)
    xCell = ds.variables['xCell']
    xEdge = ds.variables['xEdge']
    xVertex = ds.variables['xVertex']
    yCell = ds.variables['yCell']
    yEdge = ds.variables['yEdge']
    yVertex = ds.variables['yVertex']

    xMax = max(xCell)
    xMin = min(xCell)
    xMid = 0.5 * (xMin + xMax)

                                    
    yMax = max(yCell)
    yMin = min(yCell)
    yMid = 0.5 * (yMin + yMax)

    D=200.0
    R=0.02
    beta=10**(-11)
    alpha=(D/R)*beta
    A=-alpha/2.0+np.sqrt(alpha**2.0/4+(np.pi/(yMax-yMin))**2.0)
    B=-alpha/2.0-np.sqrt(alpha**2.0/4+(np.pi/(yMax-yMin))**2.0)
    Fo=1.0e-3
    g=9.8
    coriolis=2.5e-4
    gama=Fo*np.pi/(0.02*(yMax-yMin))
    p=(1.0-np.exp(B*(xMax-xMin)))/(np.exp(A*(xMax-xMin))-np.exp(B*(xMax-xMin)))#np.exp(-np.pi*xMax/yMax)
    q=1.0-p
    print(p,q,yMin)
    for iCell in range(0, nCells):
        xr=xCell[iCell]
        yr=yCell[iCell]
        
        
        velx[iCell] = gama *((yMax-yMin)/np.pi)*np.cos(np.pi*(yr-yMin)/(yMax-yMin))\
            *(p*np.exp(A*xr)+q*np.exp(B*xr)-1)

        vely[iCell] = -gama*((yMax-yMin)/np.pi)**2.0*np.sin(np.pi*(yr-yMin)/(yMax-yMin))\
                    *(p*A*np.exp(A*xr)+q*B*np.exp(B*xr))

        sshx[iCell] = -(Fo/(g*D))*((p/A)*np.exp(A*(xr-xMin))+(q/B)*np.exp(B*(xr-xMin)))-(Fo/(g*D))*((yMax-yMin)/np.pi)**2.0 \
                    *(p*A*np.exp(A*(xr-xMin))+q*B*np.exp(B*(xr-xMin)))*(np.cos(np.pi*(yr-yMin)/(yMax-yMin))-1)\
                    -((coriolis*gama/g)*((yMax-yMin)/np.pi)**2.0*np.sin(np.pi*(yr-yMin)/(yMax-yMin))-beta*(gama/g)*((yMax-yMin)/np.pi)**3.0*(np.cos(np.pi*(yr-yMin)/(yMax-yMin))-1))*(p*np.exp(A*(xr-xMin))+q*np.exp(B*(xr-xMin))-1)

    print(np.mean(sshx))
    sshx[:]=sshx[:]-np.mean(sshx[:])
    #np.savetxt('test.out',velx)

if __name__ == '__main__':
    # If called as a primary module, run main
    main()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
