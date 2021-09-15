import numpy as np

def solve(nodes, elems, forces, constraints):
    nNodes = len(nodes)
    nElems = len(elems)
    kglobal = np.zeros((3*nNodes, 3*nNodes))
    for elem in elems:
        elem.assemble(kglobal)
    kpart = kglobal
    kpart = np.delete(kpart, constraints, 0)
    kpart = np.delete(kpart, constraints, 1)
    forcesParted = np.delete(forces, constraints)
    dispVector = np.matmul(np.linalg.inv(kpart), forcesParted)
    dispTemp = np.zeros(3*nNodes)
    j = 0
    for i in range(3*nNodes):
        if i in constraints:
            pass
        else:
            dispTemp[i] = dispVector[j]
            j += 1
    dispVector = dispTemp
    forceVector = np.matmul(kglobal, dispVector)
    return dispVector, forceVector