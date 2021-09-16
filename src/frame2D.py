import numpy as np

class Frame2D:
    frame2dList = []
    def __init__(self, tag, E, A, I, node1, node2, stiffModifier = [1,1,1,1,1,1]):
        self.tag = tag
        self.E = E
        self.A = A
        self.I = I
        self.node1 = node1
        self.node2 = node2
        # a list to modify each DOF which is a list of boolean 0 or 1
        # [1,1,1,1,1,1] means beam-column element
        # [0,1,1,0,1,1] means beam element
        # you can also apply end releases with the stiffModifier
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
            self.angle = np.pi - np.arctan(dely/delx)
        elif delx < 0 and dely < 0:
            self.angle = np.pi + np.arctan(dely/delx)
        elif delx > 0 and dely < 0:
            self.angle = 2*np.pi - np.arctan(dely/delx)
        else:
            print('warning finding angle with X axix')
        Frame2D.frame2dList.append(self)



def getStiffMatrix(frame):
    A = frame.A
    E = frame.E
    L = frame.L
    I = frame.I
    ang = frame.angle
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
    # multiply stiffness modifiers to rows and columns
    for i in range(6):
        if frame.stiffModifier[i] == 0:
            ke[i, :] = 0
            ke[:, i] = 0
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

