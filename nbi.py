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
        ax.plot(np.arange(0.0,1.0,0.01),np.power(np.arange(0.0,1.0,0.01),2.0))
        ax.plot(np.arange(0.0,1.0,0.01),np.power(np.arange(0.0,1.0,0.01),3.6))
        plt.show()

def main():
    # Parse options
    op = OptionsParser()
    op = op.get_args()

    if op.type == 0:
        mps = DiscModel(op)
    elif op.type == 1:
        mps = IsoModel(op)
    elif op.type == 2:
        mps = EtaModel(op)
    elif op.type == 3:
        mps = PlummerModel(op)

    mps.create_model()

    print_model(mps.pos, mps.vel, mps.mass, op.outfile)

if __name__ == "__main__":
    main()


