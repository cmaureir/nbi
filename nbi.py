#!/usr/bin/env python

from lib.OptionsParser import *
from lib.PlummerModel  import *
from lib.DiscModel  import *
from lib.IsoModel  import *
from lib.EtaModel  import *

def plot(plot):
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(eccs,bins=40,cumulative=True,normed=True,alpha=0.5)
    #    ax.hist(eccs_small_a,bins=40,cumulative=True,normed=True,alpha=0.5)
    #    ax.hist(random_eccs,bins=10,cumulative=True,normed=True,alpha=0.5)
        ax.plot(np.arange(0.0,1.0,0.01),np.power(np.arange(0.0,1.0,0.01),2.0))
        ax.plot(np.arange(0.0,1.0,0.01),np.power(np.arange(0.0,1.0,0.01),3.6))
#        print np.arange(0.0,1.0,0.01)
    #    ax.hist(messured_a,bins=40,alpha=0.5)
        plt.show()

def main():
    # Parse options
    op = OptionsParser()
    op = op.get_args()

    disc  = op.disc
    gamma = op.gamma
    outfile = op.outfile

    #create stars
    nmp  = op.nmp
    mmp  = op.mmp
    mbh  = op.mbh
    rmax = op.ulim
    rmin = op.llim
    beta = op.beta

    if disc == 0:
        mps = DiscModel(op)
        mps.create_model()#(radii,rmin,rmax,beta,nmp,mmp,mbh)
    elif disc == 1:
        mps = IsoModel(op)
        mps.create_model()
    elif disc == 2:
        mps = EtaModel(op)
        mps.create_model()
#        mps.plot_dist()
#        mps = create_isotropic_orbits(radii,eccs,rmin,rmax,beta,nmp,mmp,mbh)
#    print mps.pos[0]
#    print mps.ID
    #    messured_eccs, messured_a = confirm_orbits(nmp,mps,mbh)
    #open file and write output

    ofile = open(outfile,'w')
    outstring = '{0:10e} {1:10e} {2:10e} {3:10e} \
                 {4:10e} {5:10e} {6:10e}\n'.format(mps.op.mbh,0.0,0.0,0.0,0.0,0.0,0.0)
    ofile.write(outstring)

    for i in range(mps.N):
        outstring = '{0:10e}\t{1:10e}\t{2:10e}\t{3:10e}\t{4:10e}\t{5:10e}\t{6:10e}\n'.format(float(mps.m[i]),\
                                                       float(mps.pos[i][0]),\
                                                       float(mps.pos[i][1]),\
                                                       float(mps.pos[i][2]),\
                                                       float(mps.vel[i][0]),\
                                                       float(mps.vel[i][1]),\
                                                       float(mps.vel[i][2]))
        ofile.write(outstring)

    ofile.close()

if __name__ == "__main__":
    main()


