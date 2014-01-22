import argparse

class OptionsParser:
    def __init__(self):

        self.parser  = None
        desc = "Create a circular, plane, stellar disc"
        self.parser = argparse.ArgumentParser(description=desc)
        #self.parser = argparse.ArgumentParser(description=desc,
        #                                      formatter_class=argparse.RawTextHelpFormatter)
        self.add_arguments()

    def add_arguments(self):
        out_help = "Name of output ascii file"
        self.parser.add_argument("-o",
                                 metavar = "outfile",
                                 dest    = "outfile",
                                 help    =  out_help,
                                 default = "disc.dat")


        ulim_help = "Outer radius for disc distribution [pc]"
        self.parser.add_argument("-u",
                                 dest     = "ulim",
                                 type     = float,
                                 help     = ulim_help,
                                 required = True)

        llim_help = "Inner radius for disc distribution [pc]"
        self.parser.add_argument("-l",
                                 dest     = "llim",
                                 type     = float,
                                 help     = llim_help,
                                 required = True)

        nmp_help = "Number of stars"
        self.parser.add_argument("-n",
                                 dest     = "nmp",
                                 type     = int,
                                 help     = nmp_help,
                                 required = True)

        mbh_help = "SMBH mass [Msun]"
        self.parser.add_argument("-mbh",
                                 dest     = "mbh",
                                 type     = float,
                                 help     = mbh_help,
                                 required = True)

        beta_help = "Surface density exponent, Sigma propto r^(-beta).\
                     If we are creating not a disk but a sphere, this should\
                     be 3D exponent - 1."
        self.parser.add_argument("-b",
                                 dest     = "beta",
                                 type     = float,
                                 help     = beta_help,
                                 required = True)

        disc_help = "0: disc,\
                     1: spherical isotropic,\
                     requires also argument for -e distribution"
        self.parser.add_argument("-d",
                                 dest    = "disc",
                                 type    = int,
                                 help    = disc_help,
                                 default = 0)

        mmp_help = "Mass of stars [Msun]"
        self.parser.add_argument("-m",
                                 dest     = "mmp",
                                 type     = float,
                                 help     = mmp_help,
                                 required = True)

        gamma_help = "Power law index gamma of eccentricity distribution,\
                      default = 1.0 (thermalized)"
        self.parser.add_argument("-e",
                                 dest     = "gamma",
                                 type     = float,
                                 help     = gamma_help,
                                 default  = 1.0,
                                 required = False)

        profile_help = "Print the density profile to this file"
        self.parser.add_argument("--profile",
                                 dest    = "profile",
                                 help    = profile_help,
                                 default = None)

        plot_help = "plot histogram of density"
        self.parser.add_argument("--plot",
                                 dest   = "plot",
                                 help   = plot_help,
                                 action = "store_true")
    def process_args(self):
        # Do something with self.parser
        pass

    def get_args(self):
        return self.parser.parse_args()
