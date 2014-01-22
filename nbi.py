#!/usr/bin/env python

from lib.OptionsParser import *
from lib.PlummerModel  import *
from lib.DiscModel  import *

def plot(plot):
    if plot:
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


def main():
    # Parse options
    op = OptionsParser()
    op = op.get_args()

    disc  = op.disc
    gamma = op.gamma
    outfile = op.outfile

    #create massive perturbers
    nmp  = op.nmp
    mmp  = op.mmp
    mbh  = op.mbh
    rmax = op.ulim
    rmin = op.llim
    beta = op.beta

    #draw semi-major axes from a power law distribution
    #the numpy.random.power uses \propto a*x^(a-1)
    radii_temp = np.random.power(2.0-beta,nmp)
    radii      = np.sort(radii_temp)

    if disc == 0:
        mps = create_star_orbits(radii,rmin,rmax,beta,nmp,mmp,mbh)

    elif disc == 1:
        eccs = np.random.power(gamma+1.0,nmp)
        mps = create_isotropic_orbits(radii,eccs,rmin,rmax,beta,nmp,mmp,mbh)

    #    messured_eccs, messured_a = confirm_orbits(nmp,mps,mbh)
    #open file and write output

    ofile = open(outfile,'w')
    outstring = '{0:10e} {1:10e} {2:10e} {3:10e} \
                 {4:10e} {5:10e} {6:10e}\n'.format(mbh,0.0,0.0,0.0,0.0,0.0,0.0)
    ofile.write(outstring)

    for mp in mps:
        outstring = '{0:10e} {1:10e} {2:10e} {3:10e} \
                     {4:10e} {5:10e} {6:10e}\n'.format(float(mp[0]),\
                                                       float(mp[1]),\
                                                       float(mp[2]),\
                                                       float(mp[3]),\
                                                       float(mp[4]),\
                                                       float(mp[5]),\
                                                       float(mp[6]))
        ofile.write(outstring)

    ofile.close()
    random_eccs_ids = random.sample(range(args.nmp),nmp)
    random_eccs = []
    for i in range(nmp):
        random_eccs.append(eccs[random_eccs_ids[i]])
    plot(op.plot, eccs)

if __name__ == "__main__":
    main()


