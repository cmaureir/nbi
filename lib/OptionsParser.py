import argparse
import textwrap

class OptionsParser:
    def __init__(self):

        # Parser
        desc = "N-body Initial conditions generator"
        self.parser = argparse.ArgumentParser(description=desc,
                                              formatter_class=argparse.RawTextHelpFormatter)
        subparsers = self.parser.add_subparsers(dest="command")

        # General options for all the models
        self.general_group = self.parser.add_argument_group("General options")
        self.add_general_arguments()

        # Convert Parser
        self.convert_parser = subparsers.add_parser("convert",
                                                 help="Convert to N-body units")
        self.add_convert_arguments()

        # Disc Parser
        self.disc_parser = subparsers.add_parser("disc",
                                                 help="Disc Model")
        self.add_disc_arguments()

        # Iso Parser
        self.iso_parser = subparsers.add_parser("iso",
                                                help="Spherical Isotropic model")
        self.add_iso_arguments()

        # Eta parser
        self.eta_parser = subparsers.add_parser("eta",
                                                help="ETA Model")
        self.add_eta_arguments()

        # Plummer parser
        self.plummer_parser = subparsers.add_parser("plummer",
                                                    help="Plummer Model")
        self.add_plummer_arguments()

    def add_convert_arguments(self):
        file_help = "Filename"
        self.convert_parser.add_argument("-f",
                                 dest     = "filename",
                                 type     = str,
                                 help     = file_help,
                                 required = False)

    def add_general_arguments(self):
        n_help = "Number of stars"
        self.general_group.add_argument("-n",
                                 dest     = "n",
                                 type     = int,
                                 help     = n_help,
                                 required = False)

        out_help = "Name of output ascii file"
        self.general_group.add_argument("-o",
                                 metavar = "outfile",
                                 dest    = "outfile",
                                 help    =  out_help,
                                 default = "nbi_output.dat")

        #profile_help = "Print the density profile to this file"
        #self.disc_parser.add_argument("--profile",
        #                         dest    = "profile",
        #                         help    = profile_help,
        #                         default = None)

        #plot_help = "plot histogram of density"
        #self.disc_parser.add_argument("--plot",
        #                         dest   = "plot",
        #                         help   = plot_help,
        #                         action = "store_true")

    def add_disc_arguments(self):

        ulim_help = "Outer radius for disc distribution [pc]"
        self.disc_parser.add_argument("-u",
                                 dest     = "ulim",
                                 type     = float,
                                 help     = ulim_help,
                                 required = True)

        llim_help = "Inner radius for disc distribution [pc]"
        self.disc_parser.add_argument("-l",
                                 dest     = "llim",
                                 type     = float,
                                 help     = llim_help,
                                 required = True)

        mbh_help = "SMBH mass [Msun]"
        self.disc_parser.add_argument("-mbh",
                                 dest     = "mbh",
                                 type     = float,
                                 help     = mbh_help,
                                 required = True)


        beta_help =  textwrap.dedent("""\
                     Surface density exponent, Sigma propto r^(-beta)""")
        self.disc_parser.add_argument("-b",
                                 dest     = "beta",
                                 type     = float,
                                 help     = beta_help,
                                 required = True)

        mmp_help = "Mass of stars [Msun]"
        self.disc_parser.add_argument("-m",
                                 dest     = "mmp",
                                 type     = float,
                                 help     = mmp_help,
                                 required = True)


    def add_eta_arguments(self):
        eta_help = "Eta parameter for eta model"
        self.eta_parser.add_argument("-eta",
                                 dest     = "eta",
                                 type     = float,
                                 help     = eta_help,
                                 required = True)

        mbh_help = "SMBH mass [Msun]"
        self.eta_parser.add_argument("-mbh",
                                 dest     = "mbh",
                                 type     = float,
                                 help     = mbh_help,
                                 required = True)

    def add_iso_arguments(self):
        ulim_help = "Outer radius for iso distribution [pc]"
        self.iso_parser.add_argument("-u",
                                 dest     = "ulim",
                                 type     = float,
                                 help     = ulim_help,
                                 required = True)


        gamma_help =  textwrap.dedent("""\
                      Power law index gamma of eccentricity distribution,
                      default = 1.0 (thermalized)""")
        self.iso_parser.add_argument("-e",
                                 dest     = "gamma",
                                 type     = float,
                                 help     = gamma_help,
                                 required = True)

        delta_help =  "3D density power law index"
        self.iso_parser.add_argument("-d",
                                 dest     = "delta",
                                 type     = float,
                                 help     = delta_help,
                                 required = True)

        mbh_help = "SMBH mass [Msun]"
        self.iso_parser.add_argument("-mbh",
                                 dest     = "mbh",
                                 type     = float,
                                 help     = mbh_help,
                                 required = True)

        mmp_help = "Mass of stars [Msun]"
        self.iso_parser.add_argument("-m",
                                 dest     = "mmp",
                                 type     = float,
                                 help     = mmp_help,
                                 required = True)


    def add_plummer_arguments(self):
        pass


    def process_args(self):
        # Do something with self.parser, maybe verify some values
        pass

    def get_args(self):
        return self.parser.parse_args()
