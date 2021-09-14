from solve import solve
import numpy as np
from frame2D import Frame2D
from node2D import Node2D
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
f3 = Frame2D(3, E, A, I, n4, n3)

deadLoad = np.zeros(3*len(Node2D.nodeList))
deadLoad[3] = -20e3
deadLoad[8] = 12e3

constraints = [0,1,2,9,10,11]


U, F = solve(Node2D.nodeList, Frame2D.frame2dList, deadLoad, constraints)
print(U)
print(F)
print('======================')


for elem in Frame2D.frame2dList:
    print('element ', elem.tag, ' displacements:')
    print(elem.getElemDisp(U))
    print('-------------------')
    print('element ', elem.tag, ' forces:')
    print(elem.getElemForce(elem.getElemDisp(U)))
    print('===================')


