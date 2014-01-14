import argparse

class OptionsParser:
    def __init__(self):
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
