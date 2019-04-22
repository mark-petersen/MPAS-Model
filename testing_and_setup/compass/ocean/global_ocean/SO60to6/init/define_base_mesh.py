#!/usr/bin/env python
"""
% Create cell width for this mesh on a regular latitude-longitude grid.
% Outputs:
%    cellWidth - m x n array, entries are desired cell width in km
%    lon - longitude, vector of length m, with entries between -180 and 180, degrees
%    lat - latitude, vector of length n, with entries between -90 and 90, degrees
"""
import numpy as np
import mesh_definition_tools as mdt


def cellWidthVsLatLon():

    ddeg = 0.1
    ddeg = 1.0

    lat = np.arange(-90, 90.01, ddeg)
    lon = np.arange(-180, 180.01, ddeg)

    cellWidthEC60to30 = mdt.EC_CellWidthVsLat(lat)
    cellWidthRRS18to6 = mdt.RRS_CellWidthVsLat(lat,18,6)

    cellWidthVsLat = mdt.mergeCellWidthVsLat(
        lat,
        cellWidthRRS18to6,
        cellWidthEC60to30,
        -60.0,  # southern transition latitude
        5.0)  # transition width

    cellWidth = np.outer(cellWidthVsLat, np.ones([1, lon.size]))

    # plot:
    #import matplotlib.pyplot as plt
    #plt.imshow(cellWidth)
    #plt.gca().invert_yaxis()
    #plt.colorbar()
    #plt.savefig('cellWidth.png')

    return cellWidth, lon, lat
