# add_initial_state.py

import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
from lxml import etree
import configparser

config = configparser.ConfigParser()
config.optionxform = lambda option: option  # preserve case for letters
config.read('config_initial_state.ini')
for section in config.sections():
    for option in config.options(section):
        try:
           vars()[option] = float(config.get(section, option))
        except:
           vars()[option] = config.get(section, option)


print(dTdx)

def ocn_generate_uniform_vertical_grid(interfaceLocations):
    nInterfaces = interfaceLocations.shape[0]
    layerSpacing = 1.0/float(nInterfaces - 1)
    interfaceLocations[0] = 0.0
    for iInterfaceLocation in range(1,nInterfaces):
        interfaceLocations[iInterfaceLocation] = interfaceLocations[iInterfaceLocation-1] + layerSpacing
    return interfaceLocations

def SpecifyInitialConditions():


    # Source file
    src = Dataset("base_mesh.nc", "r", format='NETCDF3_64BIT_OFFSET')

    # Destination file
    dst = Dataset("initial_state.nc", "w", format='NETCDF3_64BIT_OFFSET')

    # Copy attributes
    for name in src.ncattrs():
        dst.setncattr(name, src.getncattr(name))

    # Copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

    # Add dimensions
    dst.createDimension('nVertLevels', SGWnVertLevelsParam)

    srcList = list()

    # Copy variables
    for name, variable in src.variables.items():
        srcList.append(src.variables[name].name)
        x = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst.variables[name][:] = src.variables[name][:]

    # Parse the Registry file
    Registry = etree.parse("Registry.xml")
    # You can only use the above line of code after copying or linking the Registry file to this folder.
    # Alternatively, you can use the exact path location of the Registry file and modify the above line of code as
    # Registry = etree.parse("/path/to/Registry.xml")

    # Copy the units and long names of the common variables from the parsed Registry file
    for name, variable in dst.variables.items():
        for var_struct in Registry.xpath('//var_struct'):
            if var_struct.attrib['name'] == 'mesh':
                for var in var_struct.getchildren():
                    if var.attrib['name'] == name:
                        dst.variables[name].units = var.attrib['units']
                        dst.variables[name].long_name = var.attrib['description']

    # Get values of the dimensions
    nCells = len(dst.dimensions['nCells'])
    nEdges = len(dst.dimensions['nEdges'])
    nVertices = len(dst.dimensions['nVertices'])
    nVertLevels = len(dst.dimensions['nVertLevels'])

    # Create new output variables and specify their dimensions
    fEdge = dst.createVariable('fEdge', np.float64, ('nEdges',))
    fVertex = dst.createVariable('fVertex', np.float64, ('nVertices',))
    fCell = dst.createVariable('fCell', np.float64, ('nCells',))
    refZMid = dst.createVariable('refZMid', np.float64, ('nVertLevels',))
    normalVelocity = dst.createVariable('normalVelocity', np.float64, ('Time','nEdges','nVertLevels',))
    layerThickness = dst.createVariable('layerThickness', np.float64, ('Time','nCells','nVertLevels',))
    restingThickness = dst.createVariable('restingThickness', np.float64, ('nCells','nVertLevels',))
    surfaceStress = dst.createVariable('surfaceStress', np.float64, ('Time','nEdges',))
    atmosphericPressure = dst.createVariable('atmosphericPressure', np.float64, ('Time','nCells',))
    boundaryLayerDepth = dst.createVariable('boundaryLayerDepth', np.float64, ('Time','nCells',))
    refBottomDepth = dst.createVariable('refBottomDepth', np.float64, ('nVertLevels',))
    bottomDepth = dst.createVariable('bottomDepth', np.float64, ('nCells',))
    bottomDepthObserved = dst.createVariable('bottomDepthObserved', np.float64, ('nCells',))
    maxLevelCell = dst.createVariable('maxLevelCell', np.int32, ('nCells',))
    vertCoordMovementWeights = dst.createVariable('vertCoordMovementWeights', np.float64, ('nVertLevels',))
    edgeMask = dst.createVariable('edgeMask', np.int32, ('nEdges','nVertLevels',))
    temperature = dst.createVariable('temperature', np.float64, ('Time','nCells','nVertLevels',))
    salinity = dst.createVariable('salinity', np.float64, ('Time','nCells','nVertLevels',))

    # Copy the units and long names of these new output variables from the parsed Registry file
    for name, variable in dst.variables.items():
        if name not in srcList:
            for var in Registry.xpath('/registry/var_struct/var'):
                if var.attrib['name'] == name:
                    dst.variables[name].units = var.attrib['units']
                    dst.variables[name].long_name = var.attrib['description']

    # Specify the units and long names of the following variables not present in the Registry file
    temperature.units = 'degrees Celsius'
    temperature.long_name = 'potential temperature.'
    salinity.units = 'grams salt per kilogram seawater'
    salinity.long_name = 'salinity.'

    # Initialize these new output variables

    nVertLevelsP1 = nVertLevels + 1
    interfaceLocations = np.zeros(nVertLevelsP1)
    interfaceLocations = ocn_generate_uniform_vertical_grid(interfaceLocations)

    for k in range(0,nVertLevels):
        refBottomDepth[k] = SGWBottomDepthParam*interfaceLocations[k+1]
        refZMid[k] = -0.5*(interfaceLocations[k+1] + interfaceLocations[k]) \
                                       *SGWBottomDepthParam

    xCellNetCDF = src.variables['xCell']
    xCell = np.zeros(nCells)
    yCellNetCDF = src.variables['yCell']
    yCell = np.zeros(nCells)
    for iCell in range(0,nCells):
        xCell[iCell] = float(xCellNetCDF[iCell])
        yCell[iCell] = float(yCellNetCDF[iCell])

    xMin = min(xCell)
    xMax = max(xCell)
    xMidGlobal = 0.5*(xMin + xMax)
    yMin = min(yCell)
    yMax = max(yCell)
    yMidGlobal = 0.5*(yMin + yMax)

    restingThicknessArray = np.zeros((nCells,nVertLevels))
    layerThicknessArray = np.zeros((nCells,nVertLevels))

    # import from ini file later:
    constant_layer_thickness = 10.0

    for iCell in range(0,nCells):
        for k in range(0,nVertLevels):
            restingThicknessArray[iCell,k] = constant_layer_thickness
            #= SGWBottomDepthParam*(interfaceLocations[k+1] - interfaceLocations[k])
            layerThicknessArray[iCell,k] = constant_layer_thickness
            #= restingThicknessArray[iCell,k] \
            #  *(1.0 + SGWSurfaceElevationToMeanDepthRatioParam \
            #          *np.exp(-((yCell[iCell] - yMidGlobal)/SGWGaussianSurfaceElevationDecayScaleParam)**2.0))

    restingThickness[:,:] = restingThicknessArray[:,:]
    layerThickness[0,:,:] = layerThicknessArray[:,:]

    # import from ini file later:
    temperature_type = 'linear'
    T0 = 20
    dTdx = 1e-4
    dTdy = 0 #1e-4
    dTdz = 10.0/1000

    # import from ini file later:
    salinity_type = 'linear'
    S0 = 35
    dSdx = 0
    dSdy = 0
    dSdz = -10.0/1000

    for iCell in range(0,nCells):
        for k in range(0,nVertLevels):
            temperature[0,iCell,k] \
                = T0 + xCell[iCell]*dTdx + yCell[iCell]*dTdy + refZMid[k]*dTdz
            salinity[0,iCell,k] \
                = S0 + xCell[iCell]*dSdx + yCell[iCell]*dSdy + refZMid[k]*dSdz

    bottomDepth[:] = SGWBottomDepthParam
    bottomDepthObserved[:] = SGWBottomDepthParam

    fCell[:] = SGWCoriolisParam
    fEdge[:] = SGWCoriolisParam
    fVertex[:] = SGWCoriolisParam

    normalVelocity[:] = SGWNormalVelocityParam
    surfaceStress[:] = SGWSurfaceStressParam
    atmosphericPressure[:] = SGWAtmosphericPressureParam
    boundaryLayerDepth[:] = SGWBottomLayerDepthParam
    # flat bottom:
    maxLevelCell[:] = nVertLevels
    vertCoordMovementWeights[:] = SGWVertCoordMovementWeightsParam
    edgeMask[:] = SGWEdgeMaskParam

    # Close the destination file
    dst.close()

SpecifyInitialConditions()

