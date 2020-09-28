#!/usr/bin/env python
# from Sid, see https://www.dropbox.com/sh/gxo9jlvcce8ogdm/AADs16DU0MBNwrgA-_OBbbT5a/CoastalKelvinWaveMesh/ConvergenceStudyMeshes?dl=0&lst=&subfolder_nav_tracking=1

import sys
import os
import shutil
import glob
import subprocess
import numpy as np
import math

# replace with flags later
nCellsXMin = 16
nCellsXMax = 64
nCellsFactor = 2
domainWidthX = 100.0*50000.0

nCases = int(round(math.log(nCellsXMax,nCellsFactor) - math.log(nCellsXMin,nCellsFactor) + 1))
nCellsX = np.zeros(nCases,dtype=np.int8)

nCellsX[0] = nCellsXMin
nCellsFloat = nCellsXMin
for iCase in range(nCases-1):
    nCellsFloat = nCellsFloat * nCellsFactor
    # ensure nCellsX is an even integer
    nCellsX[iCase+1] = int(nCellsFloat/2+0.01)*2

for iCase in range(nCases):
    dcEdge = domainWidthX/float(nCellsX[iCase])
    os.system("planar_hex.py --nx %d --ny %d --npx --dc %f -o base_mesh_%d.nc" 
              %(nCellsX[iCase],nCellsX[iCase],dcEdge,nCellsX[iCase]))
    os.system("MpasCellCuller.x base_mesh_%d.nc culled_mesh_%d.nc" 
              %(nCellsX[iCase],nCellsX[iCase]))
    os.system("MpasMeshConverter.x culled_mesh_%d.nc mesh_%d.nc" 
              %(nCellsX[iCase],nCellsX[iCase]))
    os.system("mv culled_graph.info culled_graph_%d.info" %(nCellsX[iCase]))
    os.system("mv graph.info graph_%d.info" %(nCellsX[iCase]))

os.system("rm culled_graph_*.info")
os.system("rm graph_*.info")
