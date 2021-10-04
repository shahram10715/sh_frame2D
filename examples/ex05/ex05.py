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
n2 = Node2D(2, 3.0, 0.0)
n3 = Node2D(3, 5.0, 0.0)
n4 = Node2D(4, 7.0, 0.0)
n5 = Node2D(5, 10.0, 0.0)

E = 70e9
I = 1e-4
A = 5e-2

f1 = Frame2D(1, E, A, I, n1, n2)
f2 = Frame2D(2, E, A, I, n2, n3)
f3 = Frame2D(3, E, A, I, n3, n4)
f4 = Frame2D(4, E, A, I, n4, n5)

deadLoads = np.zeros(3*len(Node2D.nodeList))
deadLoads = applyNodeLoad(deadLoads, n3, [0.0, -80e3, 0.0])

constraints = np.zeros(3*len(Node2D.nodeList))
constraints = applyNodeConstraint(constraints, n1, [1,1,0])
constraints = applyNodeConstraint(constraints, n2, [0,1,0])
constraints = applyNodeConstraint(constraints, n4, [0,1,0])
constraints = applyNodeConstraint(constraints, n5, [0,1,0])

U, F = solve(Node2D.nodeList, Frame2D.frame2dList, deadLoads, constraints)
print(' ')
print(U)
print(F)
