#!/usr/bin/python

from lib.constants import *
from lib.GeneralRandom import *
from lib.Utils import *
import random

class Model():
    def __init__(self, n):
        self.N   = n
        self.pos = np.zeros((n,3))
        self.vel = np.zeros((n,3))
        self.mass   = np.zeros(n)
        self.ID  = np.arange(n,dtype=np.int)
