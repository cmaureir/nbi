import numpy as np
import random
import logging


def print_model(pos, vel, mass, outfile):

    # Openning file
    try:
        ofile = open(outfile,'w')
    except IOError as e:
        msg = "IO Error({0}): {1}".format(e.errno, e.strerror)
        logging.warning(msg)
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise

    # Preparing every line to print to the file
    for i in range(len(pos)):
        # Formatting particle attributes
        m  = '% 3.8e' % mass[i]
        rx = '% 3.8e' % pos[i][0]
        ry = '% 3.8e' % pos[i][1]
        rz = '% 3.8e' % pos[i][2]
        vx = '% 3.8e' % vel[i][0]
        vy = '% 3.8e' % vel[i][1]
        vz = '% 3.8e' % vel[i][2]

        # Right-align the strings
        outstring = "{0} {1} {2} {3} {4} {5} {6}\n".format( m.rjust(12),\
                                                           rx.rjust(12),\
                                                           ry.rjust(12),\
                                                           rz.rjust(12),\
                                                           vx.rjust(12),\
                                                           vy.rjust(12),\
                                                           vz.rjust(12))
        # Write to file
        ofile.write(outstring)

    # Closing the file
    ofile.close()

def spherical(r):
    """Generating 3d coordinates based on spherical coordinates"""
    vector = np.zeros(3)

    theta = np.arccos(random.uniform(-1, 1))
    phi   = random.uniform(0, 2 * np.pi)

    vector[0] = r * np.sin(theta) * np.cos(phi)
    vector[1] = r * np.sin(theta) * np.sin(phi)
    vector[2] = r * np.cos(theta)

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
                r       = pos[j] - pos[i]
                rinv    = 1 / np.linalg.norm(r)
                t_epot -= (mass[i] * mass[j]) * rinv

        t_ekin = 0.5 * mass[i] * np.dot(vel[i], vel[i])

        epot += 0.5 * t_epot
        ekin += t_ekin

    return epot, ekin

def adjust_center_of_mass(pos, vel, mass, N):
    vel_com = np.zeros(3)
    pos_com = np.zeros(3)

    for i in range(N):
        pos_com += pos[i] * mass[i]
        vel_com += vel[i] * mass[i]

    pos -= pos_com
    vel -= vel_com

def adjust_units(pos, vel, N, ekin, epot):
    alpha = -epot / 0.5
    beta  =  ekin / 0.25

    for i in range(N):
        pos[i] *= alpha
        vel[i] /= np.sqrt(beta)

def confirm_orbits(nmp,mps,mu):
    eccs  = np.zeros(nmp)
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
