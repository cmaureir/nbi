#!/usr/bin/env python

from Utils import *

class PlummerSphere:
    def __init__(self, op):

        self.mass = []
        self.pos = []
        self.vel = []

        self.epot = 0.0
        self.ekin = 0.0

        self.N = 1024
        self.M = 1.0
        self.E = -0.25
        self.r = np.array([0,0,0])
        self.v = np.array([0,0,0])
        self.a = 1

        self.scale_factor = 16.0/(3.0 * np.pi)

        self.cummulative_mass_min = 0
        self.cummulative_mass_max = 1.0/self.N

    def create_model(self):
        for i in range(self.N):
            m = self.M/self.N
            self.cummulative_mass = random.uniform(self.cummulative_mass_min, self.cummulative_mass_max)
            self.cummulative_mass_min = self.cummulative_mass_max
            self.cummulative_mass_max += m

            radius = 1 / np.sqrt(self.cummulative_mass ** (-2.0/3.0) - 1.0)

            r = spherical(radius)/self.scale_factor

            x = 0.0
            y = 0.1
            while y > (x*x) * (1.0- x*x)**3.5:
                x = random.uniform(0,1)
                y = random.uniform(0,0.1)

            velocity = x * np.sqrt(2.0) * ( 1.0 + radius * radius)**(-0.25)

            v = spherical(velocity) * np.sqrt(self.scale_factor)

            self.mass.append(m)
            self.pos.append(r)
            self.vel.append(v)

        adjust_center_of_mass(p.pos, p.vel, p.mass, p.N)
        adjust_units(p.pos, p.vel, p.N, ekin, epot)
