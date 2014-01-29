#!/usr/bin/python

from lib.constants import *
from lib.Model import *

import random
import matplotlib.pyplot as plt

class DiscModel(Model):
    def __init__(self, op):
        Model.__init__(self,op.n)

        self.mbh  = op.mbh
        self.beta = op.beta
        self.ulim = op.ulim
        self.llim = op.llim
        self.mmp  = op.mmp


    def create_model(self):
        """
        mmp = mass of perturbers
        n   = number
        mbh = black hole mass
        rmin = inner disk edge [pc]
        rmax = outer disk edge [pc]
        beta = power law index for surface density
        """
        random.seed(13223)
        # I will use the bhint format for each star
        # m x y z vx vy vz
        radii_temp = np.random.power(2.0-self.beta,self.N)
        radii      = np.sort(radii_temp)

        for i in range(self.N):
            phi = random.random() * 2.0 * np.pi
            r = radii[i] * (self.ulim-self.llim) + self.llim

            #determine keplerian velocity
            m = self.mbh
            v = np.sqrt(G*m*Msun/(r*PC_in_m)) # gives v in m/s

            #convert to km/s
            v  = v / 1000.0
            x  = np.cos(phi) * r         #*PC_in_m
            y  = np.sin(phi) * r         #*PC_in_m
            vx = -v * np.sin(phi)/9.78e5 #to convert to pc per yr
            vy =  v * np.cos(phi)/9.78e5


            self.mass[i] = self.mmp
            self.pos[i]  = np.array([x,y,0.0])
            self.vel[i]  = np.array([vx,vy,0.0])
