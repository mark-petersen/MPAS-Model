#!/usr/bin/env python
# from Sid, see https://www.dropbox.com/sh/gxo9jlvcce8ogdm/AADs16DU0MBNwrgA-_OBbbT5a/CoastalKelvinWaveMesh/ConvergenceStudyMeshes?dl=0&lst=&subfolder_nav_tracking=1

import sys
import os
import shutil
import glob
import subprocess
import numpy as np

nCellsXMin = 10
nCellsXMax = 300
d_nCellsX = 10
nCases = int((nCellsXMax - nCellsXMin)/d_nCellsX) + 1
nCellsX = np.linspace(nCellsXMin,nCellsXMax,nCases,dtype=int)

lX = 100.0*50000.0

for iCase in range(0,nCases):
    dcEdge = lX/float(nCellsX[iCase])
    os.system("./planar_hex.py --nx %d --ny %d --npx --dc %f -o base_mesh_%d.nc" 
              %(nCellsX[iCase],nCellsX[iCase],dcEdge,nCellsX[iCase]))
    os.system("MpasCellCuller.x base_mesh_%d.nc culled_mesh_%d.nc" 
              %(nCellsX[iCase],nCellsX[iCase]))
    os.system("MpasMeshConverter.x culled_mesh_%d.nc mesh_%d.nc" 
              %(nCellsX[iCase],nCellsX[iCase]))
    os.system("mv culled_graph.info culled_graph_%d.info" %(nCellsX[iCase]))
    os.system("mv graph.info graph_%d.info" %(nCellsX[iCase]))

os.system("rm culled_graph_*.info")
os.system("rm graph_*.info")
