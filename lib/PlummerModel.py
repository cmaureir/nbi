#!/usr/bin/env python

from Utils import *

class PlummerModel:
    def __init__(self, op):

        # From command line options
        self.N = op.n

        # Particle information to return after creating the model
        self.mass = []
        self.pos  = []
        self.vel  = []

        # System energy to adjust positions and velocities
        self.epot = 0.0
        self.ekin = 0.0

        # Asumptions to create Plummers spheres
        self.G = 1.0
        self.M = 1.0
        self.a = 1.0

        # Temporary position and velocity
        self.r = np.zeros(3)
        self.v = np.zeros(3)

        # Scale factor to adjust the positions and velocities
        # to have Rv = 1 (Virial radius)
        self.scale_factor = 16.0/(3.0 * np.pi)

        # Cummulative mass range to generate the radius of the particles
        self.cmass_min = 0
        self.cmass_max = 1.0/self.N

    # Function to generate the plummer sphere
    def create_model(self):

        for i in range(self.N):

            # Every particle on the system will have the same mass
            m = self.M/self.N

            # Random number between the cummulative mass range to generate
            # the radius
            self.cmass      = random.uniform(self.cmass_min, self.cmass_max)
            self.cmass_min  = self.cmass_max
            self.cmass_max += m

            # Radius constant to consider at the moment of transformation to
            # spherical coordinates
            radius = 1 / np.sqrt(self.cmass ** (-2.0/3.0) - 1.0)

            # Random position of the particle, scaled for the Rv = 1
            r = spherical(radius)/self.scale_factor

            # Random procedure to get x and y only when satisfy the
            # condition y > x^2 * (1 - x^2)^(3.5)
            x = 0.0
            y = 0.1
            while y > (x*x) * (1.0- x*x)**3.5:
                x = random.uniform(0,1)
                y = random.uniform(0,0.1)

            # Velocity constant to consider at the moment of transformation
            # to spherical coordinates
            velocity = x * np.sqrt(2.0) * ( 1.0 + radius * radius)**(-0.25)

            # Random velocity of the particle, scaled for the Rv = 1
            v = spherical(velocity) * np.sqrt(self.scale_factor)

            # Adding final information to the particle information lists
            self.mass.append(m)
            self.pos.append(r)
            self.vel.append(v)

        # Calculatin the potential and kinetic energy of the system
        self.epot, self.ekin = get_energy(self.pos, self.vel, self.mass, self.N)

        # Adjusting the center of mass of the system
        adjust_center_of_mass(self.pos, self.vel, self.mass, self.N)

        # Unit adjustment to have a total energy of -1/4 and Rv = 1
        adjust_units(self.pos, self.vel, self.N, self.ekin, self.epot)
