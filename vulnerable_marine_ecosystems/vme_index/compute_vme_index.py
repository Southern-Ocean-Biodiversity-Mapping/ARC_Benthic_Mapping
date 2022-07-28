import os
import json
import shutil
import jenkspy
import argparse
import seaborn as sns
import numpy as np
import pandas as pd
import pyreadr
from datetime import datetime
from datetime import datetime
import matplotlib.pyplot as plt
from sklearn.preprocessing import minmax_scale

# python vme_index/compute_vme_index.py -a biodata_step4.csv -v vme_index/morpho_taxa_scores_FINAL.xlsx -c vme_index/config.json -s C:/Users/cgros/code/IMAS/ARC_Data/annotation/Circumpolar_Annotation_Env_Data.RData -o 001


sns.set_style("whitegrid", {
    'grid.linestyle': '--'
 })

global LST_MORPHTAX_SMPL
LST_MORPHTAX_SMPL = [
#    "anemones_colonial-zoantharia",
#    "anemones_solitary-actiniaria",
#    "anemones_solitary-corallimorpharia",
#    "antarctic_scallop-antarctic_scallop",
    "ascidians_stalked_solitary-ascidiacea",
    "ascidians_unstalked_colonial-ascidiacea",
#    "basket_snake_stars-euryalida",
#    "brachiopods-brachiopoda",
    "bryozoans_hard_branching-bryozoans",
    "bryozoans_soft_foliaceous-bryozoans",
#    "crinoid_stalked-crinoid_stalked",
#    "hydrocorals_branching-stylasterids",
#    "hydroids_colonial_feather-hydroidolina",
#    "hydroids_solitary-hydroidolina",
    "octocorals_bottle_brush_simple-gorgonacea",
    #"octocorals_branching_bushy-alcyonacea",
    "octocorals_branching_bushy-gorgonacea",
    "octocorals_fleshy_arborescent-alcyonacea",
    "octocorals_fleshy_mushroom-alcyonacea",
    "octocorals_whip-pennatulacea",
    "octocorals_quill-pennatulacea",
#    "polychaete-serpulidae",
    "sponges_crust_encrusting-porifera",
    "sponges_erect_3d-porifera",
#    "stony_corals_solitary_free-scleractinia",
#    "stony_corals_solitary-scleractinia",
 #   "urchin_regular_pencil-cidaroida",
 #   "worms_acorn-pterobranchia"
]


def get_parser():
    parser = argparse.ArgumentParser(description='Compute VME index for each cell.',
                                     add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-a', '--afname', required=True, type=str,
                                help='Abundance data, csv file, for each cell.')
    mandatory_args.add_argument('-v', '--vfname', required=True, type=str,
                                help='Vulnerability score data, xlsx file, for each morpho-taxa.')
    mandatory_args.add_argument('-c', '--cfname', required=True, type=str,
                                help='Method parameters, JSON file.')
    mandatory_args.add_argument('-o', '--ofolder', required=True, type=str,
                                help='Output folder containing graph, config file and csv file.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-s', '--sfname', required=False, type=str,
                               help='RData filename containing sampling information.')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def quadratic_mean(x):
    vals = [xx for xx in x if str(xx) != 'nan']
    return np.sqrt(np.sum([vv * vv for vv in vals]) / len(vals))


def compute_vulnerability_score(data, agg_fct, group_avg):
    df_out = data.copy()

    if group_avg is not None:
        data["tmp"] = data[group_avg].mean(axis=1)
        lst_criteria_agg = [c for c in data.keys() if not c in ["morpho_taxon"] + group_avg]
    else:
        lst_criteria_agg = [c for c in data.keys() if not c in ["morpho_taxon"]]

    if agg_fct == "quadratic_mean":
        df_out["vulnerability_score"] = data[lst_criteria_agg].apply(quadratic_mean, axis=1)
    elif agg_fct == "sum":
        df_out["vulnerability_score"] = data[lst_criteria_agg].sum(axis=1)
    else:
        print("\nERROR: Unknown function to aggregate vulnerability scores: {} ...")
        print("\tPlease choose between quadratic_mean, and sum")
        exit()


    return df_out.sort_values(by="vulnerability_score", ignore_index=True, ascending=False)


def plot_vuln_scores(data, y_name, params, ofolder, palette, fname_out, figsize=(30,20)):
    data["taxon"] = data["morpho_taxon"].str.split('-').str[1]
    data_t_stats = data.groupby("taxon")["vulnerability_score"].aggregate(['count', 'mean', 'std', 'min', 'max']).sort_values(by="mean")
    print(data_t_stats)
    custom_dict = {}
    for idx, row in data_t_stats.iterrows():
        custom_dict[idx] = row["mean"]
    print(custom_dict)
    data = data.sort_values(by=['taxon'], key=lambda x: x.map(custom_dict))

    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(y=y_name, x="vulnerability_score", data=data,
                     orient="h", hue="taxon", ax=ax, dodge=False, palette=palette)
    title_ = ""
    for k in params.keys():
        title_ += "{} --> {}    ".format(k, params[k])
    plt.title(title_)
    plt.xlim([0, 3])

    fig.savefig(os.path.join(ofolder, fname_out), dpi=300)


def compute_abundance_score(data, lst_columns, n_breaks=0):
    df_out = data.copy()

    for col_name in lst_columns:
        if n_breaks > 0:
            breaks = jenkspy.jenks_breaks(data[data[col_name] > 0][col_name],
                                          nb_class=n_breaks)
            df_out[col_name] = pd.cut(data[col_name],
                                     bins=breaks,
                                     labels=[ll for ll in range(1, n_breaks + 1)],
                                      include_lowest=True)
            df_out[col_name].replace(to_replace=[ll for ll in range(1, n_breaks + 1)],
                                     value=[ll for ll in range(1, n_breaks + 1)], inplace=True)
            df_out.loc[data[data[col_name] == 0].index, col_name] = 0
        else:
            df_out[col_name] = minmax_scale(data[col_name], feature_range=(0, 1))

    return df_out


def plot_validation(data_abd, data_vme_idx, ofolder):
    if "VME index category" in data_vme_idx.keys():
        hue_ = data_vme_idx["VME index category"]
        palette_ = sns.color_palette("Spectral_r", max(data_vme_idx["VME index category"].unique()))
    else:
        hue_, palette_ = None, None
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.scatterplot(x=data_abd.sum(axis=1),
                    y=data_vme_idx["VME index"],
                    hue=hue_,
                    ax=ax,
                    palette=palette_)
    ax.set_xlabel("VME morpho taxa percentage cover")
    ax.set_ylabel("VME index")
    fig.savefig(os.path.join(ofolder, 'abd_v_idx.png'))
    del fig, ax

    fig, ax = plt.subplots(figsize=(10, 10))
    sns.scatterplot(x=data_abd.astype(bool).astype(int).sum(axis=1),
                    y=data_vme_idx["VME index"],
                    hue=hue_,
                    ax=ax,
                    palette=palette_)
    ax.set_xlabel("VME morpho taxa richness")
    ax.set_ylabel("VME index")
    fig.savefig(os.path.join(ofolder, 'div_v_idx.png'))
    del fig, ax


def compute_vme_index(fname_abd, fname_vuln, fname_config, folder_out, fname_sampling=None):
    print("\nLoading data ...")
    print("\tMethod parameters ...")
    with open(fname_config, 'r') as f:
        dct_config = json.load(f)
    print("\tAbundance data ...")
    df_abd = pd.read_csv(fname_abd)
    print(df_abd.head())
    print("\tVulnerability scores data ...")
    df_vuln = pd.read_excel(fname_vuln)
    print(df_vuln.head())

    folder_out = datetime.now().strftime("%Y%m%d%H%M%S") + "_" + folder_out
    print("\nCreating output folder ...")
    os.makedirs(folder_out)
    shutil.copyfile(fname_config, os.path.join(folder_out, os.path.split(fname_config)[-1]))

    lst_criteria = [c for c in df_vuln.keys() if c != "morpho_taxon"]
    df_vuln[lst_criteria] = df_vuln[lst_criteria].replace({"H": 3, "M": 2, "L": 1})
    idx_missing_scores = df_vuln[df_vuln[lst_criteria].sum(axis=1) == 0].index
    df_vuln.drop(index=idx_missing_scores, inplace=True)

    print("\nChecking if missing vulnerability scores for some morpho_taxa ...")
    lst_morpho_taxa = [mt for mt in df_abd.keys() if mt != "cellID"]
    lst_missing_scores = [mt for mt in lst_morpho_taxa if mt not in df_vuln.morpho_taxon.unique()]
    lst_missing_abd = [mt for mt in df_vuln.morpho_taxon.unique() if mt not in lst_morpho_taxa]
    df_vuln.drop(index=df_vuln[df_vuln.morpho_taxon.isin(lst_missing_abd)].index, inplace=True)
    if len(lst_missing_scores):
        print("\nMissing vulnerability scores for: {} ...".format(lst_missing_scores))
        print("\nMissing abundance data for: {} ...".format(lst_missing_abd))
        lst_morpho_taxa = [mt for mt in lst_morpho_taxa if mt not in lst_missing_scores]
    else:
        print("\tAll good ...")

    print("\nComputing the vulnerability scores using ...")
    for k in dct_config["vulnerability_score"].keys():
        print("\t{} --> {} ...".format(k, dct_config["vulnerability_score"][k]))
    df_vuln_agg = compute_vulnerability_score(data=df_vuln,
                                              agg_fct=dct_config["vulnerability_score"]["agg_fct"],
                                              group_avg=dct_config["vulnerability_score"]["group_avg"])
    print(df_vuln_agg.head(20))
    print(df_vuln_agg.tail(20))
    fname_df_vuln_score = os.path.join(folder_out, "df_vuln_score.csv")
    df_vuln_agg.to_csv(fname_df_vuln_score, index=False)
    lst_taxa = sorted(list(set([t.split("-")[-1] for t in lst_morpho_taxa])))
    palette_taxa = dict(zip(lst_taxa, sns.color_palette("Spectral", n_colors=len(lst_taxa))))
    plot_vuln_scores(df_vuln_agg,
                     y_name="morpho_taxon",
                     params=dct_config["vulnerability_score"],
                     ofolder=folder_out,
                     palette=palette_taxa,
                     fname_out="vuln_scores.png")
    plot_vuln_scores(df_vuln_agg[df_vuln_agg.morpho_taxon.isin(LST_MORPHTAX_SMPL)],
                     y_name="morpho_taxon",
                     params=dct_config["vulnerability_score"],
                     ofolder=folder_out,
                     palette=palette_taxa,
                     fname_out="vuln_scores_sample.png",
                     figsize=(15, 10))
    plot_vuln_scores(df_vuln_agg,
                     y_name="taxon",
                     params=dct_config["vulnerability_score"],
                     ofolder=folder_out,
                     palette=palette_taxa,
                     fname_out="vuln_scores_taxon.png",
                     figsize=(15, 10))

    # Abundance scores
    if dct_config["abundance_score"]["method"] == "jenks":
        print("\nCategorisation of the abundance data using Jenks natural breaks (n={}) ...".format(dct_config["abundance_score"]["param"]))
        df_abd_scr = compute_abundance_score(data=df_abd,
                                             lst_columns=lst_morpho_taxa,
                                             n_breaks=dct_config["abundance_score"]["param"])
    elif dct_config["abundance_score"]["method"] == "minmaxscale":
        print("\nScaling of the abundance data ...")
        df_abd_scr = compute_abundance_score(data=df_abd,
                                             lst_columns=lst_morpho_taxa,
                                             n_breaks=0)
    elif dct_config["abundance_score"]["method"] == "pa":
        print("\nCategorisation of the abundance data based on Presence / Absence ...")
        df_abd_scr = df_abd.copy()
        df_abd_scr[lst_morpho_taxa] = df_abd[lst_morpho_taxa].astype(bool).astype(float)
    else:
        print("\nNo categorisation of the abundance data ...")
        df_abd_scr = df_abd.copy()

    # Vulnerability index
    df_vuln_idx = df_abd_scr.copy()
    for mt in lst_morpho_taxa:
        vuln_score = df_vuln_agg[df_vuln_agg["morpho_taxon"] == mt]["vulnerability_score"].values[0]
        df_vuln_idx[mt] = df_abd_scr[mt] * vuln_score

    # VME index
    df_vme_idx = df_vuln_idx.copy()
    if dct_config["vme_index"]["agg_fct"] == "mean":
        df_vme_idx["VME index"] = df_vuln_idx[lst_morpho_taxa].mean(axis=1)
    elif dct_config["vme_index"]["agg_fct"] == "median":
        df_vme_idx["VME index"] = df_vuln_idx[lst_morpho_taxa].median(axis=1)
    elif dct_config["vme_index"]["agg_fct"] == "max":
          df_vme_idx["VME index"] = df_vuln_idx[lst_morpho_taxa].max(axis=1)
    elif dct_config["vme_index"]["agg_fct"] == "sum":
        df_vme_idx["VME index"] = df_vuln_idx[lst_morpho_taxa].sum(axis=1)
        if dct_config["abundance_score"]["method"] == "pa" and dct_config["abundance_score"]["param"]["normalize_by_area"] is True:
            df_area = pyreadr.read_r(fname_sampling)["cover_cells_env"].reset_index()[["cellID", "counts_area"]]
            df_area["cellID"] = pd.to_numeric(df_area["cellID"], downcast='integer')
            cells_no_area_info = [c for c in df_abd.cellID.to_list() if c not in df_area.cellID.to_list()]
            if len(cells_no_area_info):
                print("\tCells with no area info:")
                print(cells_no_area_info)
            df_vme_idx = pd.merge(left=df_vme_idx, right=df_area, on="cellID")
            df_vme_idx['VME index'] = df_vme_idx['VME index'] / df_vme_idx['counts_area']
    else:
        print("\nERROR: Unknown function to aggregate vulnerability indexes: {} ...")
        print("\tPlease choose between mean, median, and max")
        exit()

    # Categorisation
    if dct_config["vme_index"]["category_method"] == "jenks":
        breaks = jenkspy.jenks_breaks(df_vme_idx["VME index"], nb_class=dct_config["vme_index"]["n_breaks"])
        print(breaks)
        df_vme_idx["VME index category"] = pd.cut(df_vme_idx["VME index"],
                                                   bins=breaks,
                                                   labels=[ll for ll in range(1, dct_config["vme_index"]["n_breaks"]+1)],
                                                   include_lowest=True)
        df_vme_idx["VME index category"].replace(to_replace=[ll for ll in range(1, dct_config["vme_index"]["n_breaks"] + 1)],
                                                  value=[ll for ll in range(1, dct_config["vme_index"]["n_breaks"] + 1)],
                                                  inplace=True)
    else:
        print("WARNING: No categorisation done on the VME index ...")

    plot_validation(df_abd[lst_morpho_taxa], df_vme_idx, folder_out)

    fname_df_vme_idx = os.path.join(folder_out, "df_vme_idx.csv")
    df_vme_idx.to_csv(fname_df_vme_idx, index=False)


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    compute_vme_index(fname_abd=args.afname,
                      fname_vuln=args.vfname,
                      fname_config=args.cfname,
                      folder_out=args.ofolder,
                      fname_sampling=args.sfname)


if __name__ == "__main__":
    main()
