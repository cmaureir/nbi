#!/usr/bin/python

from lib.constants import *

import numpy as np
import sys
import argparse
import random
import matplotlib.pyplot as plt


def create_star_orbits(radii,rmin,rmax,beta,nmp,mmp,mbh):
    """
    mmp = mass of perturbers
    nmp = number
    mbh = black hole mass
    rmin = inner disk edge [pc]
    rmax = outer disk edge [pc]
    beta = power law index for surface density
    """
    random.seed(13223)
    # I will use the bhint format for each star
    # m x y z vx vy vz
    all_mps = np.zeros((nmp,7))

    for i in range(nmp):
        phi = random.random()*2.0*np.pi
        r = radii[i]*(rmax-rmin) + rmin
        #determine keplerian velocity
        mass = mbh
        v = np.sqrt(G*mass*Msun/(r*PC_in_m)) # gives v in m/s
        #convert to km/s
        v = v/1000.0
        x = np.cos(phi) * r#*PC_in_m
        y = np.sin(phi) * r#*PC_in_m
        vx = -v*np.sin(phi)/9.78e5#to convert to pc per yr
        vy = v*np.cos(phi)/9.78e5
        mpdata = np.zeros(7)
        mpdata[0] = mmp
        mpdata[1] = x
        mpdata[2] = y
        mpdata[3] = 0.0
        mpdata[4] = vx
        mpdata[5] = vy
        mpdata[6] = 0.0
        all_mps[i] = mpdata
    return all_mps

def get_pos_vel(a, b, m_anom, e, m):

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

def create_isotropic_orbits(radii,eccs,rmin,rmax,beta,nmp,mmp,mbh):
    #random.seed(13223)
    # I will use the bhint format for each star
    # m x y z vx vy vz
    all_mps = np.zeros((nmp,7))

    for i in range(nmp):
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
#        print a_vec, b_vec
        m_anomaly = random.random()*2.0*np.pi

        m, r, v = get_pos_vel(a_vec, b_vec, m_anomaly, eccs[i], mbh*Msun)
#        print m, r, v
        mpdata = np.zeros(7)
        mpdata[0] = mmp
        mpdata[1] = r[0]/PC_in_m
        mpdata[2] = r[1]/PC_in_m
        mpdata[3] = r[2]/PC_in_m
        mpdata[4] = v[0]/9.78e8
        mpdata[5] = v[1]/9.78e8
        mpdata[6] = v[2]/9.78e8
        all_mps[i] = mpdata
    return all_mps

def confirm_orbits(nmp,mps,mu):
    eccs = np.zeros(nmp)
    semis = np.zeros(nmp)
    for i in range(nmp):
        mp = mps[i]
        x = mp[1:4]
        v = mp[4:7]

#        mu = G*mu
        # r magnitude
        x_mag = np.linalg.norm(x)

        # Angular momentum
        # j = r x v
        j = np.cross(x, v)

        # Runge-Lenz-vector
        # e = { (v x j) / (G * m) }  - { r / |r| }
        e = np.cross(v, j)/(G * mu) - x/x_mag

#        print e, j, mu
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


## Main
parser = argparse.ArgumentParser(description='Create a circular, plane, stellar disc')
#parser.add_argument('-i',metavar='infile',dest='infile',help='Input file', required=True)
parser.add_argument('-o',metavar='outfile',dest='outfile',help='Name of output ascii file',default="disc.dat")
#parser.add_argument('--format',dest='format',type=int,help='Input format. 0 = Shuo, 1 = gsplash ascii. Default = 0', default=0)
parser.add_argument('-u',dest='ulim',type=float,help='Outer radius for disc distribution [pc]',required=True)
parser.add_argument('-l',dest='llim',type=float,help='Inner radius for disc distribution [pc]',required=True)
parser.add_argument('-n',dest='nmp',type=int,help='Number of stars',required=True)
parser.add_argument('-mbh',dest='mbh',type=float,help='SMBH mass [Msun]',required=True)
parser.add_argument('-b',dest='beta',type=float,help='Surface density exponent, Sigma propto r^(-beta). If we are creating not a disk but a sphere, this should be 3D exponent - 1.',required=True)
parser.add_argument('-d',dest='disc',type=int,help='0: disc,  1: spherical isotropic, requires also argument for -e distribution',default=0)
parser.add_argument('-m',dest='mmp',type=float,help='Mass of stars [Msun]',required=True)
parser.add_argument('-e',dest='gamma',type=float,help='Power law index gamma of eccentricity distribution, default = 1.0 (thermalized)',default=1.0,required=False)
parser.add_argument('--profile',dest='profile',help='Print the density profile to this file',default=None)
parser.add_argument('--plot',dest='plot',help='plot histogram of density',action='store_true')
args = parser.parse_args()

#create massive perturbers
nmp = args.nmp
mmp = args.mmp
mbh = args.mbh

rmax = args.ulim
rmin = args.llim

beta = args.beta

#draw semi-major axes from a power law distribution
#the numpy.random.power uses \propto a*x^(a-1)
radii_temp = np.random.power(2.0-beta,nmp)
radii = np.sort(radii_temp)

if args.disc == 0:
    mps = create_star_orbits(radii,rmin,rmax,beta,nmp,mmp,mbh)
elif args.disc == 1:
    eccs = np.random.power(args.gamma+1.0,nmp)
    mps = create_isotropic_orbits(radii,eccs,rmin,rmax,beta,nmp,mmp,mbh)
#    messured_eccs, messured_a = confirm_orbits(nmp,mps,mbh)
#open file and write output
ofile = open(args.outfile,'w')
outstring = '{0:10e} {1:10e} {2:10e} {3:10e} {4:10e} {5:10e} {6:10e}\n'.format(mbh,0.0,0.0,0.0,0.0,0.0,0.0)
ofile.write(outstring)
for mp in mps:
    outstring = '{0:10e} {1:10e} {2:10e} {3:10e} {4:10e} {5:10e} {6:10e}\n'.format(float(mp[0]),float(mp[1]),float(mp[2]),float(mp[3]),float(mp[4]),float(mp[5]),float(mp[6]))
    ofile.write(outstring)
ofile.close()
#eccs_small_a = []
#for i in range(len(eccs)):
#    if radii[i] < 0.00000005:
#        eccs_small_a.append(eccs[i])
random_eccs_ids = random.sample(range(args.nmp),nmp)
random_eccs = []
for i in range(nmp):
    random_eccs.append(eccs[random_eccs_ids[i]])
#print len(eccs_small_a)
if args.plot:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(eccs,bins=40,cumulative=True,normed=True,alpha=0.5)
#    ax.hist(eccs_small_a,bins=40,cumulative=True,normed=True,alpha=0.5)
#    ax.hist(random_eccs,bins=10,cumulative=True,normed=True,alpha=0.5)
    ax.plot(np.arange(0.0,1.0,0.01),np.power(np.arange(0.0,1.0,0.01),2.0))
    ax.plot(np.arange(0.0,1.0,0.01),np.power(np.arange(0.0,1.0,0.01),3.6))
    print np.arange(0.0,1.0,0.01)
#    ax.hist(messured_a,bins=40,alpha=0.5)
    plt.show()
