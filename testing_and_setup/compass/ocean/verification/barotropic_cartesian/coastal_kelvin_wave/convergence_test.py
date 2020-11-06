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
nCellsXMax = 32
nCellsFactor = math.sqrt(2)
domainWidthX = 100.0*50000.0

nCases = int(round(math.log(nCellsXMax,nCellsFactor) - math.log(nCellsXMin,nCellsFactor) + 1))
nCellsX = np.zeros(nCases,dtype=np.int8)

nCellsX[0] = nCellsXMin
nCellsFloat = nCellsXMin
for iCase in range(nCases-1):
    nCellsFloat = nCellsFloat * nCellsFactor
    # ensure nCellsX is an even integer
    nCellsX[iCase+1] = int(nCellsFloat/2+0.01)*2

print('nCellsX:',nCellsX)

for iCase in range(nCases):
    nxDir = 'nx' + str(nCellsX[iCase]).zfill(4)
    dcEdge = domainWidthX/float(nCellsX[iCase])
    os.system("mkdir %s"%(nxDir))
    os.chdir(nxDir)
    os.system("planar_hex --nx %d --ny %d --npx --dc %f -o base_mesh.nc" 
              %(nCellsX[iCase],nCellsX[iCase],dcEdge))
    os.system("MpasCellCuller.x base_mesh.nc culled_mesh.nc")
    os.system("MpasMeshConverter.x culled_mesh.nc mesh.nc")

    os.symlink('../namelist.ocean.init','namelist.ocean.init')
    os.symlink('../streams.ocean.init','streams.ocean.init')
    os.symlink('../ocean_model','ocean_model')

    # mpirun -n 4 ./ocean_model -n namelist.ocean -s streams.ocean
    subprocess.check_call(['mpirun', '-n', '1', './ocean_model', 
        '-n', 'namelist.ocean.init', '-s', 'streams.ocean.init'])

    os.symlink('../namelist.ocean.forward','namelist.ocean.forward')
    os.symlink('../streams.ocean.forward','streams.ocean.forward')
    os.symlink('../ocean_model','ocean_model')

    # mpirun -n 4 ./ocean_model -n namelist.ocean -s streams.ocean
    subprocess.check_call(['mpirun', '-n', '1', './ocean_model', 
        '-n', 'namelist.ocean.forward', '-s', 'streams.ocean.forward'])

    os.chdir("..")


