import os
import math
import copy
import argparse
import numpy as np
import pandas as pd


# Example:
#   python combine_coralnet_biigle.py -b ..\20210909_biigle_vme_cover.csv -c ..\20210909_coralnet.csv -o ..\20210915_biigle_coralnet.csv

DCT_CORALNET = {"BH_BrAnt": "bryozoans_hard_branching_antler-bryozoans",
                "BH_BrHead": "bryozoans_hard_branching_coralhead-bryozoans",
                "BH_BrLeaf": "",
                "BH_Encr": "",
                "BH_Fenestr": "",
                "BH_Lettuce": "",
                "BH_Mass": "",
                "B_Purple": "",
                "BS_Dendr": "",
                "BS_FAN": "",
                "BS_Foliac": "",
                "HydrC_Br": "",
                "S_Amorph": "",
                "S_Ball": "",
                "S_Barrel": "",
                "S_bead": "",
                "S_Buried": "",
                "S_Creep": "",
                "S_CupCmplt": "",
                "S_CupIncmp": "",
                "S_Disc": "",
                "S_Encr": "",
                "S_Er_Br": "",
                "S_Er_Lam": "",
                "S_Er_Palm": "",
                "S_Er_Simp": "",
                "S_Er_St": "",
                "S_Tube": ""
                }


def get_parser():
    parser = argparse.ArgumentParser(add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-b', '--bfname', required=True, type=str,
                                help='BIIGLE CSV filename.')
    mandatory_args.add_argument('-c', '--cfname', required=True, type=str,
                                help='CoralNet CSV filename.')
    mandatory_args.add_argument('-o', '--ofname', required=True, type=str,
                                help='CSV filename output.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def combine_coralnet_biigle(fname_biigle, fname_coralnet, fname_o):
    # Read data
    df_b = pd.read_csv(fname_biigle)
    df_c = pd.read_csv(fname_coralnet)

    # Convert CoralNet count to cover data
    tot_coralnet = df_c.sum(axis=1)
    df_c.loc[:, df_c.columns != 'Unnamed: 0'] = df_c.loc[:, df_c.columns != 'Unnamed: 0'].div(tot_coralnet, axis=0)

    # Find images in common
    print("\nTODO: FIND CORALNET METADATA and DO NOT ONLY TAKE INTERSECTION")
    list_interection_filename = list(set(df_b["filename"].tolist()) & set(df_c["Unnamed: 0"].tolist()))
    print(len(df_c), len(df_b))
    df_c = df_c[df_c["Unnamed: 0"].isin(list_interection_filename)]
    df_b = df_b[df_b["filename"].isin(list_interection_filename)]
    print(len(df_c), len(df_b))

    # Select columns of interest
    print("\nTODO: TAKE UBS_B and UBS_Sp into account?")


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    combine_coralnet_biigle(fname_biigle=args.bfname,
                            fname_coralnet=args.cfname,
                            fname_o=args.ofname)


if __name__ == "__main__":
    main()
