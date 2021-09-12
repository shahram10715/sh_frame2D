import numpy as np

class Node2D:
    nodeList = []
    def __init__(self, tag, x, y):
        self.tag = tag
        self.x = x
        self.y = y
        Node2D.nodeList.append(self)