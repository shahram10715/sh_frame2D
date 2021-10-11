import numpy as np
from node2D import *
from frame2D import *

def getGlobalMatrix(nodes, elems):
    nNodes = len(nodes)
    kglobal = np.zeros((3*nNodes, 3*nNodes))
    for elem in elems:
        kelem = getStiffMatrix(elem)
        kglobal = assemble(kglobal, elem, kelem)
    return kglobal

def solve(nodes, elems, forces, constraints):
    nNodes = len(nodes)
    fix = getFixes(constraints)
    kglobal = getGlobalMatrix(nodes, elems)
    kpart = kglobal
    kpart = np.delete(kpart, fix, 0)
    kpart = np.delete(kpart, fix, 1)
    forcesParted = np.delete(forces, fix)
    dispVector = np.matmul(np.linalg.inv(kpart), forcesParted)
    dispTemp = np.zeros(3*nNodes)
    j = 0
    for i in range(3*nNodes):
        if i not in fix:
            dispTemp[i] = dispVector[j]
            j += 1
    dispVector = dispTemp
    forceVector = np.matmul(kglobal, dispVector)
    return dispVector, forceVector

