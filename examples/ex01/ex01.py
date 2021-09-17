from solve import solve
import numpy as np
from frame2D import *
from node2D import *
###########################################################################################
#####    COPY ALL THE PY FILES IN THE MAIN DIRECTORY HERE
###########################################################################################

np.set_printoptions(precision=9)
np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=180)

n1 = Node2D(1, 0.0, 0.0)
n2 = Node2D(2, 0.0, 3.0)
n3 = Node2D(3, 4.0, 3.0)
n4 = Node2D(4, 4.0, 0.0)

E = 210e9
I = 5e-5
A = 2e-2

f1 = Frame2D(1, E, A, I, n1, n2)
f2 = Frame2D(2, E, A, I, n2, n3)
f3 = Frame2D(3, E, A, I, n3, n4)

deadLoads = np.zeros(3*len(Node2D.nodeList))
deadLoads = applyNodeLoad(deadLoads, n2, [-20e3, 0, 0])
deadLoads = applyNodeLoad(deadLoads, n3, [0, 0, 12e3])

constraints = np.zeros(3*len(Node2D.nodeList))
constraints = applyNodeConstraint(constraints, n1, [1,1,1])
constraints = applyNodeConstraint(constraints, n4, [1,1,1])

U, F = solve(Node2D.nodeList, Frame2D.frame2dList, deadLoads, constraints)

print(getElemDisp(f1, U))
print(getElemForce(f1, U))

print(getElemDisp(f2, U))
print(getElemForce(f2, U))

print(getElemDisp(f3, U))
print(getElemForce(f3, U))
