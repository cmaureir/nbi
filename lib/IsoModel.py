#!/usr/bin/python

from lib.Model import *

class IsoModel(Model):
    def __init__(self, op):
        Model.__init__(self,op.n)

        self.gamma = op.gamma
        self.delta = op.delta


    def get_pos_vel(self, a, b, m_anom, e, m):

        # Getting data
        a_vec = np.array(a)
        b_vec = np.array(b)
        m_anomaly  = m_anom
        ecc = e
        a_axis = np.linalg.norm(a_vec)
        mass = m


        # Solving the Kepler equation for elliptical orbits
        e_anomaly_new = 0
        if ecc > 0.8:
            e_anomaly_new = np.pi
        else:
            e_anomaly_new = m_anomaly

        d = 1e4
        iteration = 0
        while np.fabs(d) > 9.0e-16:
            d = e_anomaly_new - ecc * np.sin(e_anomaly_new) - m_anomaly
            if iteration - 1 >= 50:
                break
            e_anomaly_new -= d / (1.0 - ecc * np.cos(e_anomaly_new))
            iteration += 1

        e_anomaly = e_anomaly_new

        cos_e = np.cos(e_anomaly)
        sin_e = np.sin(e_anomaly)

        w = np.sqrt((G * mass)/ a_axis**3)

        r_const = cos_e - ecc
        v_const = w / (1.0 - ecc * cos_e)

        # Update position and velocity
        r = a_vec * r_const + b_vec * sin_e
        v = (-a_vec * sin_e + b_vec * cos_e) * v_const

        return mass, r, v

    def create_model(self):
        #draw semi-major axes from a power law distribution
        #the numpy.random.power uses \propto a*x^(a-1)
        eccs = np.random.power(self.gamma+1.0,self.N)

        radii_temp = np.random.power(2.0-self.op.delta,self.N)
        radii      = np.sort(radii_temp)

        for i in range(self.N):
            b = radii[i] * np.sqrt(np.fabs(1 - eccs[i]**2)) * PC_in_m
            #obtain a random 3D direction for a vector
            phi = random.random()*2.0*np.pi
            theta = random.random()*np.pi
            a_vec = np.zeros(3)
            b_vec = np.zeros(3)
            a_vec[0] = radii[i] * np.sin(theta) * np.cos(phi) * PC_in_m
            a_vec[1] = radii[i] * np.sin(theta) * np.sin(phi) * PC_in_m
            a_vec[2] = radii[i] * np.cos(theta) * PC_in_m
            #now construct a vector perpendicular to a_vec
            phi = random.random()*2.0*np.pi
            theta = random.random()*np.pi
            b_vec[0] = np.sin(theta) * np.cos(phi)
            b_vec[1] = np.sin(theta) * np.sin(phi)
            b_vec[2] = np.cos(theta) #this is now a random unit vector, right?
            #orthogonalize it to a_vec
            b_vec = b_vec - a_vec * np.dot(b_vec,a_vec) / np.dot(a_vec,a_vec)
            b_vec = b_vec/np.linalg.norm(b_vec)
            b_vec = b_vec*b
            m_anomaly = random.random()*2.0*np.pi

            m, r, v = self.get_pos_vel(a_vec, b_vec, m_anomaly, eccs[i], self.op.mbh*Msun)
            #this will write everything in PC, velocity in pc per year and mass in Msun, bhint units
            self.m[i] = self.op.mmp
            self.pos[i] = r/PC_in_m
            self.vel[i] = v[0]/9.78e8

        return
