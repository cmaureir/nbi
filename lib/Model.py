#!/usr/bin/python

import numpy as np

class Model():
    def __init__(self,op):
        self.op = op
        self.N = self.op.nmp
        self.pos = np.zeros((self.N,3))
        self.vel = np.zeros((self.N,3))
        self.m = np.zeros(self.N)
        self.ID = np.arange(self.N,dtype=np.int)
