import os
import math
import copy
import argparse
import numpy as np
import pandas as pd


# Example:
#   python ccamlr_vme_taxa_prevalence.py -i "C:/Users/cgros/data/20210806_ccamlr_records/CCAMLR_VME_Registry_09032021.xlsx" -t "C:/Users/cgros/data/20210806_ccamlr_records/taxa_scores.csv" -o ccamlr_vme_taxa_prevalence


def get_parser():
    parser = argparse.ArgumentParser(add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-i', '--ifname', required=True, type=str,
                                help='CCAMLR records excel filename.')
    mandatory_args.add_argument('-t', '--taxa', required=True, type=str,
                                help='CSV filename of the taxa vulnerability scores.')
    mandatory_args.add_argument('-o', '--ofolder', required=True, type=str,
                                help='Output folder where output are saved. If not exists, will be created.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def ccamlr_vme_taxa_prevalence(fname_i, fname_t, folder_o):
    # Read dataframe
    df = pd.read_csv(fname_i)



def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    ccamlr_vme_taxa_prevalence(fname_i=args.ifname,
                               fname_t=args.tfname,
                               folder_o=args.ofolder)


if __name__ == "__main__":
    main()
