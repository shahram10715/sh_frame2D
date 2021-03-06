import numpy as np

from node2D import Node2D

class Frame2D:
    frame2dList = []
    def __init__(self, tag, E, A, I, node1, node2, stiffModifier=[1,1]):
        self.tag = tag
        self.E = E
        self.A = A
        self.I = I
        self.node1 = node1
        self.node2 = node2
        # stiffModifier a list to apply end releases
        # [1,1] means no release
        # [1,0] means end j release
        # [0,1] means end i release
        # [0,0] means truss element
        self.stiffModifier = stiffModifier
        delx = self.node2.x - self.node1.x
        dely = self.node2.y - self.node1.y
        self.L = np.sqrt(delx**2+dely**2)
        if dely == 0 and delx > 0:
            self.angle = 0.0
        elif delx == 0 and dely > 0:
            self.angle = np.pi/2
        elif dely == 0 and delx < 0:
            self.angle = np.pi
        elif delx == 0 and dely < 0:
            self.angle = 3*np.pi/2
        elif delx > 0 and dely > 0:
            self.angle = np.arctan(dely/delx)
        elif delx < 0 and dely > 0:
            self.angle = np.pi - np.arctan(-dely/delx)
        elif delx < 0 and dely < 0:
            self.angle = np.pi + np.arctan(dely/delx)
        elif delx > 0 and dely < 0:
            self.angle = 2*np.pi - np.arctan(-dely/delx)
        else:
            print('warning finding angle with X axix')
        Frame2D.frame2dList.append(self)



def getLocalStiffMatrix(frame):
    A = frame.A
    E = frame.E
    L = frame.L
    I = frame.I
    k1 = A*E/L
    k2 = 12*E*I/L**3
    k3 = 6*E*I/L**2
    k4 = 4*E*I/L
    k5 = 2*E*I/L
    ke = np.zeros((6,6))
    ke[0, 0] = k1
    ke[3, 0] = ke[0, 3] = -k1
    ke[1, 1] = k2
    ke[2, 1] = ke[1, 2] = k3
    ke[4, 1] = ke[1, 4] = -k2
    ke[5, 1] = ke[1, 5] = k3
    ke[2, 2] = k4
    ke[4, 2] = ke[2, 4] = -k3
    ke[5, 2] = ke[2, 5] = k5
    ke[3, 3] = k1
    ke[4, 4] = k2
    ke[4, 5] = ke[5, 4] = -k3
    ke[5, 5] = k4
    # applying end releases
    if frame.stiffModifier == [1,0]:
        ind1 = [0,1,2,3,4]
        ind2 = [5]
    elif frame.stiffModifier == [0,1]:
        ind1 = [0,1,3,4,5]
        ind2 = [2]
    elif frame.stiffModifier == [0,0]:
        ind1 = [0,1,3,4]
        ind2 = [2,5]
    
    if frame.stiffModifier != [1,1]:
        k11 = []
        k22 = []
        k12 = []
        k21 = []
        for i in range(6):
            for j in range(6):
                if i in ind1 and j in ind1:
                    k11.append(ke[i,j])
                elif i in ind2 and j in ind1:
                    k21.append(ke[i,j])
                elif i in ind1 and j in ind2:
                    k12.append(ke[i,j])
                elif i in ind2 and j in ind2:
                    k22.append(ke[i,j])
        k11 = np.reshape(k11, (len(ind1), len(ind1)))
        k12 = np.reshape(k12, (len(ind1), len(ind2)))
        k21 = np.reshape(k21, (len(ind2), len(ind1)))
        k22 = np.reshape(k22, (len(ind2), len(ind2)))
        ktemp = np.linalg.inv(k22)
        ktemp = np.matmul(k12, ktemp)
        ktemp = np.matmul(ktemp, k21)
        kc = k11 - ktemp
        ke = np.zeros((6,6))
        for i in range(6):
            for j in range(6):
                if i in ind1 and j in ind1:
                    ke[i,j] = kc[ind1.index(i), ind1.index(j)]
    return ke


def getStiffMatrix(frame):
    ang = frame.angle
    ke = getLocalStiffMatrix(frame)
    R = np.zeros((6,6))
    R[0, 0] = np.cos(ang)
    R[0, 1] = np.sin(ang)
    R[1, 0] = -np.sin(ang)
    R[1, 1] = np.cos(ang)
    R[2, 2] = 1
    R[3, 3] = np.cos(ang)
    R[3, 4] = np.sin(ang)
    R[4, 3] = -np.sin(ang)
    R[4, 4] = np.cos(ang)
    R[5, 5] = 1
    ke = np.matmul((R.T), ke)
    ke = np.matmul(ke, R)
    return ke


def assemble(kglobal, frame, kelem):
    i = (frame.node1.tag-1)*3
    j = (frame.node2.tag-1)*3
    kglobal[i:i+3, i:i+3] += kelem[0:3, 0:3]
    kglobal[i:i+3, j:j+3] += kelem[0:3, 3:6]
    kglobal[j:j+3, i:i+3] += kelem[3:6, 0:3]
    kglobal[j:j+3, j:j+3] += kelem[3:6, 3:6]
    return kglobal

def getElemDisp(frame, U):
    elemDisp = np.zeros(6)
    i = 3*frame.node1.tag
    j = 3*frame.node2.tag
    elemDisp[0] = U[i-3]
    elemDisp[1] = U[i-2]
    elemDisp[2] = U[i-1]
    elemDisp[3] = U[j-3]
    elemDisp[4] = U[j-2]
    elemDisp[5] = U[j-1]
    return elemDisp

def getElemForce(frame, U):
    kelem = getStiffMatrix(frame)
    uelem = getElemDisp(frame, U)
    force = np.matmul(kelem, uelem)
    return force

