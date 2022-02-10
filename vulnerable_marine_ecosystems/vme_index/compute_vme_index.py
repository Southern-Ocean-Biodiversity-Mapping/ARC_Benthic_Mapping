import argparse
import pandas as pd

EXAMPLE_CMD = "python vme_index/compute_vme_index.py -a biodata_step4.csv -v vme_index/morpho_taxa_scores_FINAL.xlsx -o vme_index.csv"


def get_parser():
    parser = argparse.ArgumentParser(description='Compute VME index for each cell.',
                                     epilog = 'Example of use:\n\t' + EXAMPLE_CMD,
                                     add_help=True)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-a', '--afname', required=True, type=str,
                                help='Abundance data, csv file, for each cell.')
    mandatory_args.add_argument('-v', '--vfname', required=True, type=str,
                                help='Vulnerability score data, xlsx file, for each morpho-taxa.')
    mandatory_args.add_argument('-o', '--ofname', required=True, type=str,
                                help='Output data, csv file, for each cell.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def compute_vme_index(fname_abd, fname_vuln, fname_out):
    pass


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    compute_vme_index(fname_abd=args.afname,
                      fname_vuln=args.vfname,
                      fname_out=args.ofname)


if __name__ == "__main__":
    main()
