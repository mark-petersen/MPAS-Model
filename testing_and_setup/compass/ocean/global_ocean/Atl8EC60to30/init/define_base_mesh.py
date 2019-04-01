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
    cellWidthQU8 = 8 * np.ones(lat.size)
    latTransition = -30
    latWidthTransition = 5

    cellWidthAtlanticSouth = mdt.mergeCellWidthVsLat(
        lat,
        cellWidthEC60to30,
        cellWidthQU8,
        0.0,
        5.0)

    cellWidthAtlantic = mdt.mergeCellWidthVsLat(
        lat,
        cellWidthAtlanticSouth,
        cellWidthEC60to30,
        68.0,
        5.0)

    cellWidth = np.ones((lat.size, lon.size))
    for i in range(lon.size):
        cellWidth[:, i] = cellWidthVsLat

    print 'cellWidthEC60to30', cellWidthEC60to30
    print 'cellWidthQU8', cellWidthQU8
    print 'cellWidthAtlanticSouth', cellWidthAtlanticSouth
    print 'cellWidthAtlantic', cellWidthAtlantic
    #print 'cellWidthVsLat', cellWidthVsLat

    return cellWidth, lon, lat
