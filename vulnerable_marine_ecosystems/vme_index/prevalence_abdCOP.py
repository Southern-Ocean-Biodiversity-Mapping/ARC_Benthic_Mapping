import os
import shutil
import pyreadr
import argparse
import seaborn as sns
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
from compute_vme_index import LST_MORPHTAX_SMPL

# python vme_index/prevalence_abdCOP.py -i biodata_step4.csv -s biodata_sea_area.csv -o 010 -a C:/Users/cgros/code/IMAS/ARC_Data/annotation/Circumpolar_Annotation_Env_Data.RData


sns.set_style("whitegrid", {
    'grid.linestyle': '--'
 })


def get_parser():
    parser = argparse.ArgumentParser(description='Plot prevalence and abundance conditioned on presence.',
                                     add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-i', '--ifname', required=True, type=str,
                                help='Biodata, csv.')
    mandatory_args.add_argument('-s', '--sfname', required=True, type=str,
                                help='Sea area, csv.')
    mandatory_args.add_argument('-o', '--ofolder', required=True, type=str,
                                help='Output folder containing graph, config files and csv files.')
    mandatory_args.add_argument('-a', '--afname', required=True, type=str,
                                help='RData containing area / sampling effort.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def plot_prevalence(df, lst_species, folder_out, fname_out):
    df_prev_species = df[lst_species].sum(axis=0).div(len(df)).mul(100.)
    df_prev_species = df_prev_species.sort_values()
    print(df_prev_species.sort_values())

    dct_prev = {"sea": [], "species": [], "sea_color": [], "prevalence": []}
    for s in df["section_name"].unique():
        df_prev_sea = df[df["section_name"] == s]
        n = len(df_prev_sea)
        print(s)
        for t in lst_species:
            dct_prev["sea"].append(s)
            dct_prev["sea_color"].append(df_prev_sea["section_col"].unique()[0])
            dct_prev["species"].append(t)
            dct_prev["prevalence"].append(df_prev_sea[t].sum() * 100. / n)
    df_prev_plot = pd.DataFrame.from_dict(dct_prev)
    print(df_prev_plot.groupby('species')["prevalence"].std())
    print(df_prev_plot)

    fig, ax = plt.subplots(figsize=(12, 20))
    sns.barplot(data=df_prev_plot,
                order=df_prev_species.index.to_list(),
                y="species",
                x="prevalence",
                hue="sea",
                palette=dict(zip(dct_prev["sea"], dct_prev["sea_color"])),
                orient="h",
                ax=ax)
    fig.savefig(os.path.join(folder_out, fname_out), dpi=300)


def get_df_abd_plot(df, lst_species):
    dct_abd = {"species": [], "abundance_COD": []}
    for t in lst_species:
        lst_abd_t = df[df[t] > 0][t].to_list()
        dct_abd["abundance_COD"] += lst_abd_t
        dct_abd["species"] += [t for _ in lst_abd_t]
    df_abd_plot = pd.DataFrame.from_dict(dct_abd)
    return df_abd_plot


def plot_violin_abd_COD(df_abd_plot, lst_species, folder_out, suffix):
    df_abd_plot_sorted = df_abd_plot.groupby("species").max().sort_values("abundance_COD")
    print(df_abd_plot_sorted)
    lst_common_taxa = df_abd_plot_sorted.iloc[int(len(df_abd_plot_sorted) / 2) + 1:].index.to_list()
    print(lst_common_taxa)

    df_abd_t_stats = df_abd_plot.groupby('species')['abundance_COD'].aggregate(['mean', 'std', 'min', 'max']).sort_values(
        by="mean")
    print("Species abd COD")
    print(df_abd_t_stats)

    palette_taxa = dict(zip(lst_species, sns.color_palette("Spectral", n_colors=len(lst_species))))

    fig, ax = plt.subplots(figsize=(13, 5))
    sns.violinplot(data=df_abd_plot[df_abd_plot["species"].isin(lst_common_taxa)],
                   x="species",
                   y="abundance_COD",
                   scale="count",
                   order=df_abd_t_stats.loc[lst_common_taxa].index.to_list(),
                   palette=palette_taxa,
                   orient="v",
                   cut=0,
                   inner=None,
                   bw="scott",
                   ax=ax)
    ax.set_ylim(-0.2, 30)
    #plt.xticks(rotation=45)
    fig.savefig(os.path.join(folder_out, 'abd_COP_common'+suffix+'.png'), dpi=300)

    fig, ax = plt.subplots(figsize=(13, 5))
    sns.violinplot(data=df_abd_plot[~df_abd_plot["species"].isin(lst_common_taxa)],
                   x="species",
                   y="abundance_COD",
                   scale="count",
                   order=df_abd_t_stats.loc[
                       [t for t in df_abd_t_stats.index.to_list() if t not in lst_common_taxa]].index.to_list(),
                   palette=palette_taxa,
                   orient="v",
                   cut=0,
                   inner=None,
                   bw="scott",
                   ax=ax)
    ax.set_ylim(-0.02, 4)
    #plt.xticks(rotation=45)
    fig.savefig(os.path.join(folder_out, 'abd_COP_rare_'+suffix+'.png'), dpi=300)


def prevalence_abdCOP(fname_input, fname_sea, folder_out, fname_sampling):
    print("\nLoading data ...")
    print("\tInput data ...")
    df = pd.read_csv(fname_input)
    print("\tDiversity data ...")
    df_sea = pd.read_csv(fname_sea)
    print(df_sea.section_name.unique())

    lst_morpho_taxa = [mt for mt in df.keys() if mt != "cellID"]
    #lst_morpho_taxa = LST_MORPHTAX_SMPL
    lst_taxa = sorted(list(set([t.split("-")[-1] for t in lst_morpho_taxa])))
    for t in lst_taxa:
        lst_morpho_cur = [mt for mt in lst_morpho_taxa if mt.endswith(t)]
        df[t] = df[lst_morpho_cur].sum(axis=1)

    folder_out = datetime.now().strftime("%Y%m%d%H%M%S") + "_" + folder_out
    print("\nCreating output folder ...")
    os.makedirs(folder_out)

    # Bar plot prevalence
    df_prev = df[["cellID"] + lst_taxa + lst_morpho_taxa]
    df_prev[lst_taxa + lst_morpho_taxa] = df_prev[lst_taxa + lst_morpho_taxa].astype(bool).astype(int)
    df_prev = pd.merge(df_prev, df_sea, on="cellID", how="left")

    #plot_prevalence(df=df_prev,
    #                lst_species=LST_MORPHTAX_SMPL,
    #                folder_out=folder_out,
    #                fname_out="bar_plot_smpl_mt.png")

    # Violin plot diversity

    # Divide the richness by the sampling effort of the cell
    df_area = pyreadr.read_r(fname_sampling)["cover_cells_env"].reset_index()[["cellID", "counts_area"]]
    df_area["cellID"] = pd.to_numeric(df_area["cellID"], downcast='integer')
    df_prev = pd.merge(left=df_prev, right=df_area, on="cellID")

    dct_div = {"sea": [], "morphotaxa_richness": [], "sea_color": []}
    for s in df_prev["section_name"].unique():
        df_prev_sea = df_prev[df_prev["section_name"] == s]
        df_prev_sea_pa = df_prev_sea[lst_morpho_taxa].astype(bool).astype(float)
        df_prev_sea_pa_count = df_prev_sea_pa.sum(axis=1) / df_prev_sea["counts_area"]
        dct_div["sea"] += df_prev_sea["section_name"].to_list()
        dct_div["sea_color"] += df_prev_sea["section_col"].to_list()
        dct_div["morphotaxa_richness"] += df_prev_sea_pa_count.to_list()

    df_div_plot = pd.DataFrame.from_dict(dct_div)
    print("Morpho taxa richness mean {:.2f} std {:.2f} max {} ...".format(df_div_plot["morphotaxa_richness"].mean(),
                                                                          df_div_plot["morphotaxa_richness"].std(),
                                                                          df_div_plot["morphotaxa_richness"].max()))

    df_div_mt_stats = df_div_plot.groupby('sea')['morphotaxa_richness'].aggregate(['mean', 'std', 'min', 'max']).sort_values(by="mean")
    print("Morpho taxa richness")
    print(df_div_mt_stats)

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.violinplot(data=df_div_plot,
                   x="morphotaxa_richness",
                   y="sea",
                   scale="count",
                   order=df_div_mt_stats.index.to_list(),
                   palette=dict(zip(dct_div["sea"], dct_div["sea_color"])),
                   orient="h",
                   cut=0,
                   inner=None,
                   bw="scott",
                   ax=ax)
    ax.set_xlim(-0.2, 15)
    fig.savefig(os.path.join(folder_out, 'violin_mt_div.png'), dpi=300)
    exit()
    dct_div = {"sea": [], "taxa_richness": [], "sea_color": []}
    for s in df_prev["section_name"].unique():
        df_prev_sea = df_prev[df_prev["section_name"] == s]
        df_prev_sea_pa = df_prev_sea[lst_taxa].astype(bool).astype(float)
        df_prev_sea_pa_count = df_prev_sea_pa.sum(axis=1)
        dct_div["sea"] += df_prev_sea["section_name"].to_list()
        dct_div["sea_color"] += df_prev_sea["section_col"].to_list()
        dct_div["taxa_richness"] += df_prev_sea_pa_count.to_list()

    df_div_plot = pd.DataFrame.from_dict(dct_div)

    print("Taxa richness mean {:.2f} std {:.2f} max {} ...".format(df_div_plot["taxa_richness"].mean(),
                                                                   df_div_plot["taxa_richness"].std(),
                                                                   df_div_plot["taxa_richness"].max()))

    df_div_t_stats = df_div_plot.groupby('sea')['taxa_richness'].aggregate(['mean', 'std', 'min', 'max']).sort_values(by="mean")
    print("Taxa richness")
    print(df_div_t_stats)

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.violinplot(data=df_div_plot,
                   x="taxa_richness",
                   y="sea",
                   scale="count",
                   order=df_div_t_stats.index.to_list(),
                   palette=dict(zip(dct_div["sea"], dct_div["sea_color"])),
                   orient="h",
                   cut=0,
                   inner=None,
                   bw="scott",
                   ax=ax)
    fig.savefig(os.path.join(folder_out, 'violin_t_div.png'), dpi=300)

    # Violin plot abundance COD
    df_abd_plot_taxa = get_df_abd_plot(df, lst_taxa)
    df_abd_plot_mt = get_df_abd_plot(df, LST_MORPHTAX_SMPL)

    plot_violin_abd_COD(df_abd_plot_taxa, lst_taxa, folder_out, "taxa")
    plot_violin_abd_COD(df_abd_plot_mt, LST_MORPHTAX_SMPL, folder_out, "mt")

    df_stats_taxa = df_abd_plot_taxa.groupby('species')['abundance_COD'].aggregate(['count', 'mean', 'std', 'min', 'max'])
    df_stats_taxa["prevalence"] = df_stats_taxa["count"].div(len(df)).mul(100.)
    df_stats_taxa.drop(columns="count", inplace=True)
    df_stats_taxa.to_csv(os.path.join(folder_out, "stats_taxa.csv"), float_format='%.3f')

    dct_abd = {"morpho_taxa": [], "abundance_COD": []}
    for mt in lst_morpho_taxa:
        lst_abd_mt = df[df[mt] > 0][mt].to_list()
        dct_abd["abundance_COD"] += lst_abd_mt
        dct_abd["morpho_taxa"] += [mt for _ in lst_abd_mt]
    df_abd_plot = pd.DataFrame.from_dict(dct_abd)

    df_stats_morphotaxa = df_abd_plot.groupby('morpho_taxa')['abundance_COD'].aggregate(['count', 'mean', 'std', 'min', 'max'])
    df_stats_morphotaxa["prevalence"] = df_stats_morphotaxa["count"].div(len(df)).mul(100.)
    df_stats_morphotaxa.drop(columns="count", inplace=True)
    df_stats_morphotaxa.to_csv(os.path.join(folder_out, "stats_morphotaxa.csv"), float_format='%.3f')

    shutil.copyfile(fname_input, os.path.join(folder_out, "biodata_step4.csv"))


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    prevalence_abdCOP(fname_input=args.ifname,
                      fname_sea=args.sfname,
                      folder_out=args.ofolder,
                      fname_sampling=args.afname)


if __name__ == "__main__":
    main()
