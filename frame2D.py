import numpy as np
from node2D import Node2D


class Frame2D:
    frame2dList = []
    def __init__(self, E, A, I, node1, node2):
        self.E = E
        self.A = A
        self.I = I
        self.node1 = node1
        self.node2 = node2
        self.L = np.sqrt((self.node1.x-self.node2.x)**2+(self.node1.y-self.node2.y)**2)
        delx = self.node2.x - self.node1.x
        dely = self.node2.y - self.node1.y
        vec1 = np.array([delx, dely])
        vec1 = vec1/np.linalg.norm(vec1)
        vec2 = np.array([1,0])
        dp = np.dot(vec1, vec2)
        self.angle = np.arccos(dp)
        Frame2D.frame2dList.append(self)

    def getStiffMatrix(self):
        A = self.A
        E = self.E
        L = self.L
        I = self.I
        ang = self.angle
        k1 = A*E/L
        k2 = 12*E*I/L**3
        k3 = 6*E*I/L**2
        k4 = 4*E*I/L
        k5 = 2*E*I/L
        ke = np.zeros((6,6))
        ke[0][0] = k1
        ke[3][0] = ke[0][3] = -k1
        ke[1][1] = k2
        ke[2][1] = ke[1][2] = k3
        ke[4][1] = ke[1][4] = -k2
        ke[5][1] = ke[1][5] = k3
        ke[2][2] = k4
        ke[4][2] = ke[2][4] = -k3
        ke[5][2] = ke[2][5] = k5
        ke[3][3] = k1
        ke[4][4] = k2
        ke[4][5] = ke[5][4] = -k3
        ke[5][5] = k4
        R = np.zeros((6,6))
        R[0][0] = np.cos(ang)
        R[0][1] = np.sin(ang)
        R[1][0] = -np.sin(ang)
        R[1][1] = np.cos(ang)
        R[2][2] = 1
        R[3][3] = np.cos(ang)
        R[3][4] = np.sin(ang)
        R[4][3] = -np.sin(ang)
        R[4][4] = np.cos(ang)
        R[5][5] = 1
        ke = np.matmul((R.T), ke)
        ke = np.matmul(ke, R)
        return ke

    def assemble(self, kglobal):
        i = self.node1.tag
        j = self.node2.tag
        kelem = self.getStiffMatrix()
        kglobal[3*i-3][3*i-3] += kelem[0][0]
        kglobal[3*i-3][3*i-2] += kelem[0][1]
        kglobal[3*i-3][3*i-1] += kelem[0][2]
        kglobal[3*i-3][3*j-3] += kelem[0][3]
        kglobal[3*i-3][3*j-2] += kelem[0][4]
        kglobal[3*i-3][3*j-1] += kelem[0][5]

        kglobal[3*i-2][3*i-3] += kelem[1][0]
        kglobal[3*i-2][3*i-2] += kelem[1][1]
        kglobal[3*i-2][3*i-1] += kelem[1][2]
        kglobal[3*i-2][3*j-3] += kelem[1][3]
        kglobal[3*i-2][3*j-2] += kelem[1][4]
        kglobal[3*i-2][3*j-1] += kelem[1][5]

        kglobal[3*i-1][3*i-3] += kelem[2][0]
        kglobal[3*i-1][3*i-2] += kelem[2][1]
        kglobal[3*i-1][3*i-1] += kelem[2][2]
        kglobal[3*i-1][3*j-3] += kelem[2][3]
        kglobal[3*i-1][3*j-2] += kelem[2][4]
        kglobal[3*i-1][3*j-1] += kelem[2][5]

        kglobal[3*j-3][3*i-3] += kelem[3][0]
        kglobal[3*j-3][3*i-2] += kelem[3][1]
        kglobal[3*j-3][3*i-1] += kelem[3][2]
        kglobal[3*j-3][3*j-3] += kelem[3][3]
        kglobal[3*j-3][3*j-2] += kelem[3][4]
        kglobal[3*j-3][3*j-1] += kelem[3][5]

        kglobal[3*j-2][3*i-3] += kelem[4][0]
        kglobal[3*j-2][3*i-2] += kelem[4][1]
        kglobal[3*j-2][3*i-1] += kelem[4][2]
        kglobal[3*j-2][3*j-3] += kelem[4][3]
        kglobal[3*j-2][3*j-2] += kelem[4][4]
        kglobal[3*j-2][3*j-1] += kelem[4][5]

        kglobal[3*j-1][3*i-3] += kelem[5][0]
        kglobal[3*j-1][3*i-2] += kelem[5][1]
        kglobal[3*j-1][3*i-1] += kelem[5][2]
        kglobal[3*j-1][3*j-3] += kelem[5][3]
        kglobal[3*j-1][3*j-2] += kelem[5][4]
        kglobal[3*j-1][3*j-1] += kelem[5][5]

        return kglobal


# test main problem here
np.set_printoptions(precision=4)
n1 = Node2D(1, 0.0, 0.0)
n2 = Node2D(2, 0.0, 3.0)
n3 = Node2D(3, 4.0, 3.0)
n4 = Node2D(4, 4.0, 0.0)

E = 210e9
I = 5e-5
A = 2e-2

f1 = Frame2D(E, A, I, n1, n2)
f2 = Frame2D(E, A, I, n2, n3)
f3 = Frame2D(E, A, I, n4, n3)


kglobal = np.zeros((12,12))
#breakpoint()
f1.assemble(kglobal)
f2.assemble(kglobal)
f3.assemble(kglobal)
print(kglobal[5])



