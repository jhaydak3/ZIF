#!/usr/bin/env python 
import numpy as np
import sys
from readFiles import *

thisfh = sys.argv[1]
linkerfh = "part_Im.xyz"

#Read the linker file
lAtomList, lAtomCord = readxyz(linkerfh)
sAtomList, sAtomCord = readxyz(thisfh)
a,b,c,alpha,beta,gamma = readcifFile(thisfh[:-4] + ".cif")
cell_params = [a, b, c, alpha, beta, gamma]

#sAtomCord, sAtomList = reduceToUnitCell(sAtomCord,sAtomList,cell_params,-1,2)

sAtomList = replaceSiwZn(sAtomList)
#writexyzFile(sAtomCord,sAtomList,"testZn.xyz")

minDist = calcMinZnZnDist(sAtomCord,sAtomList)
sf = 6/minDist
a = a*sf; b = b*sf; c = c*sf;

sAtomCord = expandStructure(sAtomCord,sf)
#writexyzFile(sAtomCord,sAtomList,"testZnExpanded.xyz")

sAtomCord, sAtomList = putLinkerIn(sAtomList,lAtomList,sAtomCord,lAtomCord)

cell_params = [a, b, c, alpha, beta, gamma]
#writexyzFile(sAtomCord,sAtomList, thisfh[:-4] + "_ZIF.xyz",cell_params)	

reducedCord,reducedList = reduceToUnitCell(sAtomCord,sAtomList,cell_params,.5,1.5)
writexyzFile(reducedCord,reducedList, thisfh[:-4] + "_ZIF_unitcell.xyz",cell_params)
