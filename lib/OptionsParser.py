import argparse
import textwrap

class OptionsParser:
    def __init__(self):

        self.parser  = None
        desc = "N-body Initial conditions generator"
        #desc = "Create a circular, plane, stellar disc"
        #self.parser = argparse.ArgumentParser(description=desc)
        self.parser = argparse.ArgumentParser(description=desc,
                                              formatter_class=argparse.RawTextHelpFormatter)

        self.general_group = self.parser.add_argument_group("General options")
        self.add_general_arguments()

        self.disc_group = self.parser.add_argument_group("Disc model")
        self.add_disc_arguments()

        self.iso_group = self.parser.add_argument_group("Spherical Isotropic model")
        self.add_iso_arguments()

        self.eta_group = self.parser.add_argument_group("ETA model")
        self.add_eta_arguments()

        self.plummer_group = self.parser.add_argument_group("Plummer model")
        self.add_plummer_arguments()

    def add_general_arguments(self):
        model_help =  textwrap.dedent("""\
                       0: Disc,
                       1: Spherical Isotropic,
                       2: Eta model,
                       3: Plummer Model""")
        self.general_group.add_argument("-t",
                                 dest    = "type",
                                 type    = int,
                                 help    = model_help,
                                 default = 0,
                                 required = True)

        n_help = "Number of stars"
        self.general_group.add_argument("-n",
                                 dest     = "n",
                                 type     = int,
                                 help     = n_help,
                                 required = True)

        #out_help = "Name of output ascii file"
        #self.general_group.add_argument("-o",
        #                         metavar = "outfile",
        #                         dest    = "outfile",
        #                         help    =  out_help,
        #                         default = "disc.dat")

        #profile_help = "Print the density profile to this file"
        #self.disc_group.add_argument("--profile",
        #                         dest    = "profile",
        #                         help    = profile_help,
        #                         default = None)

        #plot_help = "plot histogram of density"
        #self.disc_group.add_argument("--plot",
        #                         dest   = "plot",
        #                         help   = plot_help,
        #                         action = "store_true")

    def add_disc_arguments(self):

        ulim_help = "Outer radius for disc distribution [pc]"
        self.disc_group.add_argument("-u",
                                 dest     = "ulim",
                                 type     = float,
                                 help     = ulim_help,
                                 required = False)

        llim_help = "Inner radius for disc distribution [pc]"
        self.disc_group.add_argument("-l",
                                 dest     = "llim",
                                 type     = float,
                                 help     = llim_help,
                                 required = False)

        mbh_help = "SMBH mass [Msun]"
        self.disc_group.add_argument("-mbh",
                                 dest     = "mbh",
                                 type     = float,
                                 help     = mbh_help,
                                 required = False)


        beta_help =  textwrap.dedent("""\
                     Surface density exponent, Sigma propto r^(-beta)""")
        self.disc_group.add_argument("-b",
                                 dest     = "beta",
                                 type     = float,
                                 help     = beta_help,
                                 required = False)

        mmp_help = "Mass of stars [Msun]"
        self.disc_group.add_argument("-m",
                                 dest     = "mmp",
                                 type     = float,
                                 help     = mmp_help,
                                 required = False)


    def add_eta_arguments(self):
        eta_help = "Eta parameter for eta model"
        self.eta_group.add_argument("-eta",
                                 dest     = "eta",
                                 type     = float,
                                 help     = eta_help,
                                 required = False)

    def add_iso_arguments(self):
        gamma_help =  textwrap.dedent("""\
                      Power law index gamma of eccentricity distribution,
                      default = 1.0 (thermalized)""")
        self.disc_group.add_argument("-e",
                                 dest     = "gamma",
                                 type     = float,
                                 help     = gamma_help,
                                 default  = 1.0,
                                 required = False)

    def add_plummer_arguments(self):
        pass


    def process_args(self):
        # Do something with self.parser, maybe verify some values
        pass

    def get_args(self):
        return self.parser.parse_args()
