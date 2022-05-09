import os
import math
import copy
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# Example:
#   python ccamlr_vme_taxa_prevalence.py -i "C:/Users/cgros/data/20210806_ccamlr_records/CCAMLR_VME_Registry_09032021.xlsx" -t "C:/Users/cgros/data/20210806_ccamlr_records/taxa_scores.csv" -o ccamlr_vme_taxa_prevalence


def get_parser():
    parser = argparse.ArgumentParser(add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-i', '--ifname', required=True, type=str,
                                help='CCAMLR records excel filename.')
    mandatory_args.add_argument('-t', '--tfname', required=True, type=str,
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
    xls = pd.ExcelFile(fname_i)
    df = pd.read_excel(xls, 'VME risk areas (taxa)', skiprows=2)

    if os.path.isdir(folder_o):
        print("Output folder already exists {} ...".format(folder_o))
    else:
        print("Creating output folder {} ...".format(folder_o))
        os.makedirs(folder_o)

    # Split name columns
    df["taxon_name"] = df["VMEIndicatorTaxon"].str.split('(').str[0]
    df["taxon_code"] = df["VMEIndicatorTaxon"].str.split('(').str[1].str.split(')').str[0]

    print("\nTODO: Weight plots\n")
    metric = "VMESpecimenCount"
    #metric = "VMESpecimenWeight"

    df_sum = df[["taxon_code", metric]].groupby(['taxon_code']).sum().reset_index()

    df_taxa = pd.read_csv(fname_t)
    list_taxa = df_taxa["Taxon_Code"].tolist()

    dct_not_recorded = {"taxon_code": [], metric: []}
    for t in list_taxa:
        if t not in df_sum["taxon_code"].tolist():
            dct_not_recorded["taxon_code"].append(t)
            dct_not_recorded[metric].append(0)
    df_not_recorded = pd.DataFrame.from_dict(dct_not_recorded)
    df_sum = pd.concat([df_sum, df_not_recorded])
    df_sum.sort_values(metric, inplace=True)

    palette = sns.color_palette("hls", len(df_sum))

    fname_o_total = os.path.join(folder_o, metric+"_tot.png")
    plt.figure(figsize=(20, 10))
    sns.barplot(data=df_sum, x="taxon_code", y=metric, palette=palette)
    plt.title("Total {} for each VME taxon.".format(metric))
    plt.xlabel("Taxon code")
    plt.ylabel(metric)
    print("Saving plot {} ...".format(fname_o_total))
    plt.savefig(fname_o_total)
    #plt.show()

    fname_o_median = os.path.join(folder_o, metric + "_median.png")
    plt.figure(figsize=(20, 10))
    sns.barplot(data=df, x="taxon_code", y=metric, estimator=np.median, ci=95, palette=palette, order=df_sum["taxon_code"].tolist())
    plt.title("Median {} for each VME taxon across Risk Areas (95% CI).".format(metric))
    plt.xlabel("Taxon code")
    plt.ylabel(metric)
    print("Saving plot {} ...".format(fname_o_median))
    plt.savefig(fname_o_median)
    #plt.show()



def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    ccamlr_vme_taxa_prevalence(fname_i=args.ifname,
                               fname_t=args.tfname,
                               folder_o=args.ofolder)


if __name__ == "__main__":
    main()
