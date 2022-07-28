import os
import shutil
import pyreadr
import argparse
import seaborn as sns
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt


# python vme_index/correlation_area_richness.py -d 20220301164405_005 -a C:\Users\cgros\code\IMAS\ARC_Data\annotation\Circumpolar_Annotation_Env_Data.RData -o 005


sns.set_style("whitegrid", {
    'grid.linestyle': '--'
 })


def get_parser():
    parser = argparse.ArgumentParser(description='Assess the relation between Sampling effort and Diversity-based VME indexes.',
                                     add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-d', '--dfname', required=True, type=str,
                                help='Diversity VME index data, folder with csv and config file.')
    mandatory_args.add_argument('-a', '--afname', required=True, type=str,
                                help='Rdata metadata file.')
    mandatory_args.add_argument('-o', '--ofolder', required=True, type=str,
                                help='Output folder containing graph, config files and csv files.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def correlation_area_richness(folder_richness, fname_area, folder_out):
    print("\nLoading data ...")
    print("\tArea data ...")
    df_area = pyreadr.read_r(fname_area)["cell_metadata_env"].reset_index()
    df_area = df_area[["cellID", "counts_area"]]
    print("\tDiversity data ...")
    fname_rich = os.path.join(folder_richness, "df_vme_idx.csv")
    df_rich = pd.read_csv(fname_rich)
    df_rich = df_rich[["cellID", "VME index"]]
    df_rich.rename(columns={"VME index": "Diversity VME index"}, inplace=True)

    df_ = pd.merge(df_area, df_rich, on="cellID", how="right")
    print(df_.head(10))

    folder_out = datetime.now().strftime("%Y%m%d%H%M%S") + "_" + folder_out
    print("\nCreating output folder ...")
    os.makedirs(folder_out)

    fig, ax = plt.subplots(figsize=(8, 8))
    sns.scatterplot(data=df_, x="counts_area", y="Diversity VME index")
    fig.savefig(os.path.join(folder_out, 'corr_area_richness.png'), dpi=300)

    shutil.copyfile(fname_rich, os.path.join(folder_out, "df_vme_idx_div.csv"))
    shutil.copyfile(os.path.join(folder_richness, "config.json"), os.path.join(folder_out, "config_div.json"))


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    correlation_area_richness(folder_richness=args.dfname,
                                   fname_area=args.afname,
                                   folder_out=args.ofolder)


if __name__ == "__main__":
    main()
