import os
import math
import copy
import argparse
import numpy as np
import pandas as pd
import pyreadr


# Example:
#   python data_preparation\combine_coralnet_biigle.py -b biodata_step2.csv -c biodata_step3.csv -m C:\Users\cgros\code\IMAS\ARC_Data\annotation\Circumpolar_Annotation_Data.Rdata -a coverage_biigle839 -o 20211116


LST_BIIGLE839_FULL = ["PS81_shallow", "PS61", "PS14", "PS06", "JR17001", "JR262", "AA2011", "JR15005", "PS96", "JR17003"]
LST_PS81SHALLOW_TRANSECT = ["185", "186", "189"]


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
    mandatory_args.add_argument('-a', '--afolder', required=True, type=str,
                                help='Folder containing the area covered by BIIGLE 839 annotations.')
    mandatory_args.add_argument('-o', '--ofolder', required=True, type=str,
                                help='Output folder.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def combine_coralnet_biigle(fname_biigle, fname_coralnet, fname_metadata, folder_area_notfull, folder_o):
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
    df_m.loc[:, "biigle839"] = [0 for i in df_m.coralnet.tolist()]
    df_m.drop(index=df_m[df_m.coralnet == 0].index, inplace=True)

    # Fill BIIGLE839
    # Get infos
    #exit()
    df_m["survey"] = df_m["filename"].str.split('_').str[0]
    df_m["transect"] = df_m["filename"].str.split('_').str[1]
    df_m["imageID"] = df_m["filename"].str.split('_').str[2]
    # Fill 1s the surveys that have been fully annotated
    for survey_full in LST_BIIGLE839_FULL:
        if survey_full != "PS81_shallow":
            idx_survey_full = df_m[df_m["survey"] == survey_full].index
            df_m.loc[idx_survey_full, "biigle839"] = 1
        else:
            for transect_ps81shallow in LST_PS81SHALLOW_TRANSECT:
                idx_transect_full = df_m[(df_m["survey"] == "PS81") & (df_m["transect"] == transect_ps81shallow)].index
                df_m.loc[idx_transect_full, "biigle839"] = 1
    # Fill 1s the surveys that have been partially annotated
    lst_biigle839_notfull = [s for s in list(set(df_m["survey"].tolist())) if s not in LST_BIIGLE839_FULL]
    for survey_notfull in lst_biigle839_notfull:
        fname_survey_notfull = os.path.join(folder_area_notfull, survey_notfull+".xlsx")
        if os.path.isfile(fname_survey_notfull):
            df_survey_notfull = pd.read_excel(fname_survey_notfull, dtype=str)
            for transect_notfull in list(set(df_survey_notfull["Transect"].tolist())):
                idx_transect_notfull = df_m[(df_m["survey"] == survey_notfull) & (df_m["transect"] == transect_notfull)].index
                if len(idx_transect_notfull):
                    start = df_survey_notfull[df_survey_notfull["Transect"] == transect_notfull]["Start"].tolist()[0]
                    end = df_survey_notfull[df_survey_notfull["Transect"] == transect_notfull]["End"].tolist()[0]
                    start_int, end_int = int(start), int(end)
                    list_annotated_str = [str(z).zfill(len(start)) for z in list(range(start_int, end_int+1))]
                    idx_annotated = df_m[(df_m["survey"] == survey_notfull)
                                       & (df_m["transect"] == transect_notfull)
                                       & (df_m["imageID"].isin(list_annotated_str))].index
                    df_m.loc[idx_annotated, "biigle839"] = 1
                elif (survey_notfull == "PS81" and transect_notfull == "159") \
                        or (survey_notfull == "CRS" and transect_notfull == "1103")\
                        or (survey_notfull == "tan1901" and transect_notfull == "209")\
                        or (survey_notfull == "TAN1802" and transect_notfull in ["180", "160", "195", "94", "213",
                                                                                 "196", "184", "191", "208", "179",
                                                                                 "207", "209", "183", "92", "97",
                                                                                 "193", "170", "98", "197", "96",
                                                                                 "185"]):
                    pass
                    # Ignoring these transects because of bad image quality
                else:
                    print(survey_notfull, transect_notfull)
                    # Issue here
                    exit()
        else:
            print("ERROR: file not found: {} ...".format(fname_survey_notfull))
            exit()

    df_m.drop(columns=["survey", "transect", "imageID"], inplace=True)
    print("\nNumber of images annotated with:")
    print("\tCoralNet: {} ...".format(len(df_m[df_m["coralnet"] == 1].index)))
    print("\tBIIGLE254: {} ...".format(len(df_m[df_m["biigle254"] == 1].index)))
    print("\tBIIGLE839: {} ...".format(len(df_m[df_m["biigle839"] == 1].index)))

    list_area_m = list(set([row["filename"] for i_r, row in df_m.iterrows() if row["area"] != "nan"]))
    list_area_b = list(set([row["filename"] for i_r, row in df_b.iterrows() if row["area"] != "nan"]))
    if len([f for f in list_area_b if f not in list_area_m]):
        print("ERROR: BIIGLE has more area data than METADATA.")
        print([f for f in list_area_b if f not in list_area_m][:10])

    print("\nTODO: Check missing Area missing data.")
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

    df_m["survey"] = df_m["filename"].str.split('_').str[0]

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
                            folder_area_notfull=args.afolder,
                            folder_o=args.ofolder)


if __name__ == "__main__":
    main()
