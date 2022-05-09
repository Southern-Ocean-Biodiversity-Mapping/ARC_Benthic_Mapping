import os
import json
import shutil
import jenkspy
import argparse
import seaborn as sns
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
from sklearn.preprocessing import minmax_scale

# python vme_index/correlation_abundance_richness.py -a 20220218162214_002 -d 20220218161751_002 -o 004


sns.set_style("whitegrid", {
    'grid.linestyle': '--'
 })


def get_parser():
    parser = argparse.ArgumentParser(description='Assess the relation between the Abundance- and Diversity-based VME indexes.',
                                     add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-a', '--afname', required=True, type=str,
                                help='Abundance VME index data, folder with csv and config file.')
    mandatory_args.add_argument('-d', '--dfname', required=True, type=str,
                                help='Diversity VME index data, folder with csv and config file.')
    mandatory_args.add_argument('-o', '--ofolder', required=True, type=str,
                                help='Output folder containing graph, config files and csv files.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def correlation_abundance_richness(folder_abd, folder_richness, folder_out):
    print("\nLoading data ...")
    print("\tAbundance data ...")
    fname_abd = os.path.join(folder_abd, "df_vme_idx.csv")
    df_abd = pd.read_csv(fname_abd)
    df_abd = df_abd[["cellID", "VME index"]]
    df_abd.rename(columns={"VME index": "Abundance VME index"}, inplace=True)
    print("\tDiversity data ...")
    fname_rich = os.path.join(folder_richness, "df_vme_idx.csv")
    df_rich = pd.read_csv(fname_rich)
    df_rich = df_rich[["cellID", "VME index"]]
    df_rich.rename(columns={"VME index": "Diversity VME index"}, inplace=True)

    df_ = pd.merge(df_abd, df_rich, on="cellID", how="right")
    print(df_.head(10))

    folder_out = datetime.now().strftime("%Y%m%d%H%M%S") + "_" + folder_out
    print("\nCreating output folder ...")
    os.makedirs(folder_out)

    fig, ax = plt.subplots(figsize=(10, 10))
    sns.scatterplot(data=df_, x="Abundance VME index", y="Diversity VME index")
    fig.savefig(os.path.join(folder_out, 'assemblage_count.png'), dpi=300)

    shutil.copyfile(fname_abd, os.path.join(folder_out, "df_vme_idx_abd.csv"))
    shutil.copyfile(fname_rich, os.path.join(folder_out, "df_vme_idx_div.csv"))
    shutil.copyfile(os.path.join(folder_abd, "config.json"), os.path.join(folder_out, "config_abd.json"))
    shutil.copyfile(os.path.join(folder_richness, "config.json"), os.path.join(folder_out, "config_div.json"))


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    correlation_abundance_richness(folder_abd=args.afname,
                                   folder_richness=args.dfname,
                                   folder_out=args.ofolder)


if __name__ == "__main__":
    main()
