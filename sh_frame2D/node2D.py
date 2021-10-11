import numpy as np

class Node2D:
    nodeList = []
    def __init__(self, tag, x, y):
        self.tag = tag
        self.x = x
        self.y = y
        Node2D.nodeList.append(self)


def applyNodeLoad(loads, node, nodeLoad):
    loads[3*node.tag-3] = nodeLoad[0]
    loads[3*node.tag-2] = nodeLoad[1]
    loads[3*node.tag-1] = nodeLoad[2]
    return loads

def applyNodeConstraint(constraints, node, nodeConstraint):
    constraints[3*node.tag-3] = nodeConstraint[0]
    constraints[3*node.tag-2] = nodeConstraint[1]
    constraints[3*node.tag-1] = nodeConstraint[2]
    return constraints

def getFixes(constraints):
    fix = []
    for i in range(len(constraints)):
        if constraints[i] == 1:
            fix.append(i)
    return fix


