#!/usr/bin/python

from lib.constants import *
from lib.Model import *
from lib.GeneralRandom import *
from lib.Utils import *

import numpy as np
import random

import matplotlib.pyplot as plt

class EtaModel(Model):
    def __init__(self, op):
        Model.__init__(self,op.n)

        self.mbh  = op.mbh
        self.eta = op.eta

    def distribution_function(self,r):
        return self.eta/(4.0*np.pi)/(np.power(r,3.0-self.eta) * np.power((1.0+r),self.eta-1.0))

    def create_model(self):
        #draw semi-major axes from a power law distribution
        #the numpy.random.power uses \propto a*x^(a-1)
        mbh = self.mbh
        #total mass in stars should be unity
        mmp = 1.0/(self.N) #Equal mass stars
        self.eta = self.eta
        self.rtemp = np.arange(0.01,1.0,0.001)
        self.p = self.distribution_function(self.rtemp)
        dist = GeneralRandom(x=self.rtemp,p=self.p)

#        fig = plt.figure()
#        ax = fig.add_subplot(111)
#        ax.plot(dist.x,dist.pdf,label='pdf')
#        ax.plot(dist.x,dist.cdf,label='cdf')
#        ax.plot(dist.y,dist.inversecdf,label='inv cdf')#,normed=True,bins=50)
#        ax.legend()
#        plt.show()
#

        self.radii = dist.random(self.N)[0]
        #now I have the radii, obtain vel dispersions for each radius
        self.veldisp = np.power(self.radii,3.0-self.eta)*np.power((1.0+self.radii),1.0+self.eta)\
            *(1.0/(2.0*self.eta-4.0)-4.0/(2.0*self.eta-3.0)+3.0/(self.eta-1.0)-4.0/(2.0*self.eta-1.0)+0.5*self.eta)\
            -np.power(self.radii,self.eta-1.0)*np.power(self.radii+1.0,5.0-self.eta)/(2.0*self.eta-4.0)\
            +np.power(self.radii,self.eta)*4.0*np.power(self.radii+1.0,4.0-self.eta)/(2.0*self.eta-3.0)\
            -3.0*np.power(self.radii,self.eta+1.0)*np.power(self.radii+1.0,3.0-self.eta)/(self.eta-1.0)\
            +4.0*np.power(self.radii,self.eta+2.0)*np.power(self.radii+1.0,2.0-self.eta)/(2.0*self.eta-1.0)\
            -np.power(self.radii,self.eta+3.0)*np.power(self.radii+1.0,1.0-self.eta)/(2.0*self.eta)
        #for each star, determine three independent gaussian random variables with this sigma std
        for i in range(self.N):
            self.vel[i] = np.random.normal(0.0,self.veldisp[i],size=3)
            self.pos[i] = spherical(self.radii[i])
            self.mass[i] = mmp

        return

    def plot_dist(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.rtemp,self.p)
        ax.hist(self.radii[0],normed=True,bins=50)
        plt.show()

        return
