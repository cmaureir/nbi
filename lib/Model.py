#!/usr/bin/python

import numpy as np

class Model():
    def __init__(self, n):
        self.N   = n
        self.pos = np.zeros((n,3))
        self.vel = np.zeros((n,3))
        self.mass   = np.zeros(n)
        self.ID  = np.arange(n,dtype=np.int)
