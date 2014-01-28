#!/usr/bin/env python

import random
import numpy as np


class PlummerSphere:
    def __init__(self):

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

    def spherical(self, r):

        vector = np.array([0.0, 0.0, 0.0])
        theta = np.arccos(random.uniform(-1, 1))
        phi = random.uniform(0, 2*np.pi)

        vector[0] = r * np.sin( theta ) * np.cos( phi )
        vector[1] = r * np.sin( theta ) * np.sin( phi )
        vector[2] = r * np.cos( theta )

        return vector


    def create_model(self):
        for i in range(self.N):
            m = self.M/self.N
            self.cummulative_mass = random.uniform(self.cummulative_mass_min, self.cummulative_mass_max)
            self.cummulative_mass_min = self.cummulative_mass_max
            self.cummulative_mass_max += m

            radius = 1 / np.sqrt(self.cummulative_mass ** (-2.0/3.0) - 1.0)

            r = self.spherical(radius)/self.scale_factor

            x = 0.0
            y = 0.1
            while y > (x*x) * (1.0- x*x)**3.5:
                x = random.uniform(0,1)
                y = random.uniform(0,0.1)

            velocity = x * np.sqrt(2.0) * ( 1.0 + radius * radius)**(-0.25)

            v = self.spherical(velocity) * np.sqrt(self.scale_factor)

            #print m, r[0], r[1], r[2], v[0], v[1], v[2]
            self.mass.append(m)
            self.pos.append(r)
            self.vel.append(v)





if __name__ == "__main__":
    p = PlummerSphere()
    p.create_model()
    p.get_energy(p.pos, p.vel,p.mass)
    p.adjust_center_of_mass()
    print "--------------"
    p.get_energy(p.pos, p.vel,p.mass)
    p.adjust_units()
    print "--------------"
    p.get_energy(p.pos, p.vel,p.mass)
