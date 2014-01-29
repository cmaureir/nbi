import numpy as np
import random


def print_model(pos, vel, mass, outfile):

    ofile = open(outfile,'w')
    #outstring = '{0:10e} {1:10e} {2:10e} {3:10e} \
    #             {4:10e} {5:10e} {6:10e}\n'.format(mbh,0.0,0.0,0.0,0.0,0.0,0.0)
    #ofile.write(outstring)

    for i in range(len(pos)):
        outstring = '{0:10e}\t{1:10e}\t{2:10e}\t{3:10e}\t{4:10e}\t{5:10e}\t{6:10e}\n'.format(float(mass[i]),\
                                                       float(pos[i][0]),\
                                                       float(pos[i][1]),\
                                                       float(pos[i][2]),\
                                                       float(vel[i][0]),\
                                                       float(vel[i][1]),\
                                                       float(vel[i][2]))
        ofile.write(outstring)

    ofile.close()

def spherical(r):

    vector = np.array([0.0, 0.0, 0.0])
    theta = np.arccos(random.uniform(-1, 1))
    phi = random.uniform(0, 2*np.pi)

    vector[0] = r * np.sin( theta ) * np.cos( phi )
    vector[1] = r * np.sin( theta ) * np.sin( phi )
    vector[2] = r * np.cos( theta )

    return vector


def get_energy(pos, vel, mass, N):

    # Total Potential and Kinetic Energy
    epot = 0.0
    ekin = 0.0

    # Temporal values
    t_epot = 0.0
    t_ekin = 0.0

    for i in range(N):
        t_epot = 0
        for j in range(N):
            if i != j:
                rx = pos[j][0] - pos[i][0]
                ry = pos[j][1] - pos[i][1]
                rz = pos[j][2] - pos[i][2]
                rinv = 1 / np.sqrt(rx*rx + ry*ry + rz*rz)

                t_epot -= (mass[i] * mass[j]) * rinv
        t_ekin = 0.5 * mass[i] * (vel[i][0]**2 + vel[i][1]**2 + vel[i][2]**2)

        epot += 0.5 * t_epot
        ekin += t_ekin

    return epot, ekin

def adjust_center_of_mass(pos, vel, mass, N):
    vel_com = np.zeros(3)
    pos_com = np.zeros(3)

    for i in range(N):
        pos_com[0] += pos[i][0] * mass[i]
        pos_com[1] += pos[i][1] * mass[i]
        pos_com[2] += pos[i][2] * mass[i]

        vel_com[0] += vel[i][0] * mass[i]
        vel_com[1] += vel[i][1] * mass[i]
        vel_com[2] += vel[i][2] * mass[i]

    for i in range(N):
        pos[i][0] -= pos_com[0]
        pos[i][1] -= pos_com[1]
        pos[i][2] -= pos_com[2]

        vel[i][0] -= vel_com[0]
        vel[i][1] -= vel_com[1]
        vel[i][2] -= vel_com[2]


def adjust_units(pos, vel, N, ekin, epot):
    alpha = -epot / 0.5
    beta  =  ekin / 0.25

    for i in range(N):
        pos[i] *= alpha
        vel[i] /= np.sqrt(beta)

def confirm_orbits(nmp,mps,mu):
    eccs = np.zeros(nmp)
    semis = np.zeros(nmp)
    for i in range(nmp):
        mp = mps[i]
        x = mp[1:4]
        v = mp[4:7]

        # r magnitude
        x_mag = np.linalg.norm(x)

        # Angular momentum
        # j = r x v
        j = np.cross(x, v)

        # Runge-Lenz-vector
        # e = { (v x j) / (G * m) }  - { r / |r| }
        e = np.cross(v, j)/(G * mu) - x/x_mag

        # Eccentricity
        ecc = np.linalg.norm(e)

        # Semi-major axis
        # a = ( j * j ) / (G * m * | 1 - ecc^2 | )
        a = j.dot(j)/(G* mu * np.fabs(1-ecc**2))

        semis[i] = a
        eccs[i] = ecc
    return eccs, semis

def print_profile(m,r,fname):
    """Print the density profile

    Arguments:
    - `m`:
    - `r`:
    - `fname`:
    """
    #find out what dr is
    ofile = open(fname,'w')
    dr = r[1] - r[0]
    volfac = dr*4.0*np.pi
    for i in range(1,len(r)):
        radius = r[i]
        mass = m[i]
        dens = mass /( radius*radius*volfac)
        outline = "{0:.6e} {1:.6e}\n".format(radius,dens)
        ofile.write(outline)
    ofile.close()
