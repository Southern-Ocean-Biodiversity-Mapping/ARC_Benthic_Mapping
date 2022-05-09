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

# python vme_index/assemblage_uniqueness.py -a biodata_step4.csv -o 003


sns.set_style("whitegrid", {
    'grid.linestyle': '--'
 })


def get_parser():
    parser = argparse.ArgumentParser(description='Quantify the uniqueness of VME indicator taxa assemblages.',
                                     add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-a', '--afname', required=True, type=str,
                                help='Abundance data, csv file, for each cell.')
    mandatory_args.add_argument('-o', '--ofolder', required=True, type=str,
                                help='Output folder containing graph, config file and csv file.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def plot_count_per_cell(data, ofolder, palette, title):
    fig, ax = plt.subplots(figsize=(35, 20))
    sns.countplot(y="assemblage", data=data, palette=palette)
    ax.set_xlabel("Count of cells per assemblage")
    plt.title(title)
    fig.savefig(os.path.join(ofolder, 'assemblage_count.png'), dpi=300)


def compute_assemblage_uniqueness(fname_abd, folder_out):
    print("\nLoading data ...")
    print("\tAbundance data ...")
    df_abd = pd.read_csv(fname_abd)
    print(df_abd.head())

    folder_out = datetime.now().strftime("%Y%m%d%H%M%S") + "_" + folder_out
    print("\nCreating output folder ...")
    os.makedirs(folder_out)

    lst_taxa = sorted(list(set([t.split("-")[1] for t in df_abd.keys() if t != "cellID"])))
    for taxa_ in lst_taxa:
        lst_morpho_taxa_ = list(set([t for t in df_abd.keys() if t.endswith(taxa_)]))
        df_abd[taxa_] = df_abd[lst_morpho_taxa_].sum(axis=1)

    print("\nCategorisation of the abundance data based on Presence / Absence ...")
    df_abd_pa = df_abd.copy()[["cellID"] + lst_taxa]
    df_abd_pa[lst_taxa] = df_abd[lst_taxa].astype(bool).astype(float)

    print("\nGet assemblages for each cell ...")
    for idx, row in df_abd_pa.iterrows():
        lst_taxa_cell = sorted(list(set([c for c in lst_taxa if row[c] == 1])))
        df_abd_pa.loc[idx, "assemblage"] = "__".join(lst_taxa_cell)

    print(df_abd_pa.head())
    lst_assemblage = sorted(list(set(df_abd_pa.assemblage.unique())))

    dct_assemblage = {"assemblage": [], "count": []}
    for a in lst_assemblage:
        dct_assemblage["assemblage"].append(a)
        dct_assemblage["count"].append(len(df_abd_pa[df_abd_pa.assemblage == a]))
    df_assemblage = pd.DataFrame.from_dict(dct_assemblage)
    df_assemblage.sort_values(by="count", inplace=True)

    n_assemblage = len(df_assemblage) - 1
    print("\nNumber of different taxonomic assemblages: {} ...".format(n_assemblage))
    n_cell_unique = len(df_assemblage[df_assemblage["count"] == 1])
    print("\nNumber of cells where a given assemblage occurs only in this cell: {} ... "
          "\n\t{:.2f} % of cells) ... "
          "\n\t{:.2f} % of the assemblages occur only in one cell ...".format(n_cell_unique, n_cell_unique * 100. / len(df_abd_pa), n_cell_unique * 100. / n_assemblage))
    n_cell_empty = df_assemblage[df_assemblage["assemblage"] == ""]["count"].values[0]
    print("\nNumber of cells where there is no VME indicator taxon: {} ({:.2f} % of cells) ...".format(n_cell_empty, n_cell_empty * 100. / len(df_abd_pa)))

    lst_common_assemblage = df_assemblage[df_assemblage["count"] >= 5]["assemblage"].to_list()
    lst_common_assemblage = [a for a in lst_common_assemblage if a != ""]
    palette_common_assemblage = dict(zip(lst_common_assemblage, sns.color_palette("Spectral", n_colors=len(lst_common_assemblage))))
    plot_count_per_cell(data=df_abd_pa[df_abd_pa.assemblage.isin(lst_common_assemblage)],
                        ofolder=folder_out,
                        palette=palette_common_assemblage,
                        title="Assemblage occuring in 5 cells or more")


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    compute_assemblage_uniqueness(fname_abd=args.afname, folder_out=args.ofolder)


if __name__ == "__main__":
    main()
