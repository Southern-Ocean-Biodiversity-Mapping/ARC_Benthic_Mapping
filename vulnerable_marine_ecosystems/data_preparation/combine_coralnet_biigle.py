import os
import math
import copy
import argparse
import numpy as np
import pandas as pd
import pyreadr


# Example:
#   python data_preparation\combine_coralnet_biigle.py -b 20210917_biigle_vme_cover.csv -c 20210917_coralnet_cover.csv -m ../../ARC_data/Circumpolar_Annotation_Data.Rdata -o 20210920


def get_parser():
    parser = argparse.ArgumentParser(add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-b', '--bfname', required=True, type=str,
                                help='BIIGLE CSV filename.')
    mandatory_args.add_argument('-c', '--cfname', required=True, type=str,
                                help='CoralNet CSV filename.')
    mandatory_args.add_argument('-m', '--mfname', required=True, type=str,
                                help='Metadata Rdata filename.')
    mandatory_args.add_argument('-o', '--ofolder', required=True, type=str,
                                help='Output folder.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def combine_coralnet_biigle(fname_biigle, fname_coralnet, fname_metadata, folder_o):
    # Read data
    df_b = pd.read_csv(fname_biigle)
    df_c = pd.read_csv(fname_coralnet)
    df_m = pyreadr.read_r(fname_metadata)["image_metadata"].reset_index()

    # Cleanup
    df_m.drop(columns=["rownames", "proj_coord_x", "proj_coord_y", "cellID"], inplace=True)
    df_m.rename(columns={"Filename.standardised": "filename",
                         "lon": "longitude",
                         "lat": "latitude",
                         "cover": "coralnet",
                         "counts": "biigle254"}, inplace=True)
    df_m.loc[:, "coralnet"] = [1 if i == "yes" else 0 for i in df_m.coralnet.tolist()]
    df_m.loc[:, "biigle254"] = [1 if i == "yes" else 0 for i in df_m.biigle254.tolist()]
    df_m.loc[:, "biigle839"] = [0 for i in df_m.biigle254.tolist()]
    df_m.drop(index=df_m[df_m.coralnet == 0].index, inplace=True)
    print("\nTODO: Pull biigle839 annotated images.")

    list_area_m = list(set([row["filename"] for i_r, row in df_m.iterrows() if row["area"] != "nan"]))
    list_area_b = list(set([row["filename"] for i_r, row in df_b.iterrows() if row["area"] != "nan"]))
    if len([f for f in list_area_b if f not in list_area_m]):
        print("ERROR: BIIGLE has more area data than METADATA.")
        print([f for f in list_area_b if f not in list_area_m][:10])

    print(df_b.keys())
    print(df_c.keys())
    print(df_m.keys())

    print("\nTODO: Check missing Area missing data + disprecancy between BIIGLE vs Jan's data.")
    df_b.drop(columns=["area"], inplace=True)

    dct_taxon_source = {"morpho_taxon": [], "source": []}
    for taxon_coralnet in [c for c in df_c.keys() if c not in ["filename", "n_annotation"]]:
        dct_taxon_source["morpho_taxon"].append(taxon_coralnet)
        dct_taxon_source["source"].append("coralnet")
    for taxon_biigle839 in [c for c in df_b.keys() if c not in ['filename', 'area', 'area_pix']]:
        dct_taxon_source["morpho_taxon"].append(taxon_biigle839)
        dct_taxon_source["source"].append("biigle839")
    #for taxon_biigle254 in [c for c in df_b.keys() if c not in ['filename', 'area', 'area_pix']]:
    #    dct_taxon_source["morpho_taxon"].append(taxon_biigle254)
    #         dct_taxon_source["source"].append("biigle254")
    df_taxon_source = pd.DataFrame.from_dict(dct_taxon_source)
    df_m.drop(columns=["coralnet", "biigle839", "biigle254"], inplace=True)
    print("\nTODO: Pull Jan's and Victor's data.")

    if len(df_m[df_m["longitude"].isnull()]):
        print("ERROR: MISSING LONGITUDE")
        # df_nan_longitude_rows = df[df["image_longitude"].isnull()]
        # print("Dropping {} rows with missing longitude info...".format(len(df_nan_longitude_rows)))
        # print("\tTODO: FETCH THIS INFO...")
        # df.drop(df_nan_longitude_rows.index, axis="index", inplace=True)
        exit()
    if len(df_m[df_m["latitude"].isnull()]):
        print("ERROR: MISSING LATITUDE")
        exit()

    print("\nMerging datasets ...")
    df_merged = pd.merge(df_m, df_b, on="filename", how="left")
    df_merged = pd.merge(df_merged, df_c, on="filename", how="left")

    if not os.path.isdir(folder_o):
        print("\nCreating output folder: {} ...".format(folder_o))
        os.makedirs(folder_o)

    fname_o = os.path.join(folder_o, "bio_data.csv")
    print("\nSaving data: {} ...".format(fname_o))
    df_merged.to_csv(fname_o, index=False)
    fname_o_src = os.path.join(folder_o, "bio_data_source.csv")
    print("\nSaving data: {} ...".format(fname_o_src))
    df_taxon_source.to_csv(fname_o_src, index=False)


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    combine_coralnet_biigle(fname_biigle=args.bfname,
                            fname_coralnet=args.cfname,
                            fname_metadata=args.mfname,
                            folder_o=args.ofolder)


if __name__ == "__main__":
    main()
