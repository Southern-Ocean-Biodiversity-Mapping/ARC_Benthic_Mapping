import os
import argparse
import pandas as pd
import pyreadr

import curate_biigle_report


# Example:
#   python data_preparation\combine_coralnet_biigle.py -b biodata_step2.csv -c biodata_step3.csv -m C:\Users\cgros\code\IMAS\ARC_Data\annotation\Circumpolar_Annotation_Data.Rdata -a coverage_biigle839 -p px_size.csv -o biodata_step4.csv


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
    mandatory_args.add_argument('-p', '--pfname', required=True, type=str,
                                help='File containing the area, in sqpix, of each image.')
    mandatory_args.add_argument('-o', '--ofname', required=True, type=str,
                                help='Output csv filename.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def combine_coralnet_biigle(fname_biigle, fname_coralnet, fname_metadata, folder_area_notfull, fname_sqpx, fname_o):
    # Read data
    df_b = pd.read_csv(fname_biigle)
    df_c = pd.read_csv(fname_coralnet)
    df_m = pyreadr.read_r(fname_metadata)["image_metadata"].reset_index()
    df_p = pd.read_csv(fname_sqpx)[["filename", "image_size_sqpx"]]

    print("\nWARNING: If add AREA: carefull of summing differently between CoralNet BIIGLE839 BIIGLE254.\n")

    # Cleanup
    df_m.drop(columns=["rownames", "proj_coord_x", "proj_coord_y", 'area_source', 'image_quality_score', 'year', 'gear', 'area'], inplace=True)
    df_m.rename(columns={"Filename.standardised": "filename",
                         "lon": "longitude",
                         "lat": "latitude",
                         "cover": "coralnet",
                         "counts": "biigle254"}, inplace=True)
    df_m.loc[:, "coralnet"] = [1 if i == "yes" else 0 for i in df_m.coralnet.tolist()]
    df_m.loc[:, "biigle254"] = [1 if i == "yes" else 0 for i in df_m.biigle254.tolist()]
    df_m.loc[:, "biigle839"] = [0 for i in df_m.coralnet.tolist()]
    df_m.drop(index=df_m[df_m.coralnet == 0].index, inplace=True)
    df_m["cellID"] = df_m["cellID"].astype(int)

    # Fill BIIGLE839
    # Get infos
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

    df_m.drop(columns=["transect", "imageID"], inplace=True)

    #print("\nRemoving {} images because area is unknown ...".format(len(df_m[df_m.area.isnull()])))
    #filename_to_exclude = df_m[df_m.area.isnull()]["filename"].tolist()
    #df_m = df_m[~df_m.filename.isin(filename_to_exclude)]
    #df_b = df_b[~df_b.filename.isin(filename_to_exclude)]
    #df_c = df_c[~df_c.filename.isin(filename_to_exclude)]

    filename_to_exclude = df_c[df_c.n_annotation == 0]["filename"].tolist()
    if len(filename_to_exclude):
        print("\nRemoving {} images because no annotation on CoralNet ...".format(len(df_c[df_c.n_annotation == 0])))
        df_m = df_m[~df_m.filename.isin(filename_to_exclude)]
        df_b = df_b[~df_b.filename.isin(filename_to_exclude)]
        df_c = df_c[~df_c.filename.isin(filename_to_exclude)]

    print("\nNumber of images annotated with:")
    print("\tCoralNet: {} ...".format(len(df_m[df_m["coralnet"] == 1].index)))
    print("\tBIIGLE254: {} ...".format(len(df_m[df_m["biigle254"] == 1].index)))
    print("\tBIIGLE839: {} ...".format(len(df_m[df_m["biigle839"] == 1].index)))

    list_im_biigle = list(set(df_m[df_m.biigle254 == 1]["filename"].tolist() + df_m[df_m.biigle839 == 1]["filename"].tolist()))
    list_ann_biigle = list(set(df_b["filename"].tolist()))
    list_no_ann_biigle = [f for f in list_im_biigle if f not in list_ann_biigle]
    print("\nAdding zeros to biigle report where no annotation was made, {} images ...".format(len(list_no_ann_biigle)))
    dct_no_ann = {"filename": list_no_ann_biigle}
    for taxon_biigle in [c for c in df_b.keys() if c not in ['filename']]:
        dct_no_ann[taxon_biigle] = [0 for _ in list_no_ann_biigle]
    df_b = pd.concat([df_b, pd.DataFrame.from_dict(dct_no_ann)])

    df_b = pd.merge(df_p, df_b, on="filename", how="right")

    df_c = pd.merge(df_m[['cellID', 'filename']], df_c, on="filename", how="right")
    df_b = pd.merge(df_m[['cellID', 'filename']], df_b, on="filename", how="right")

    lst_taxa_coralnet = [t for t in df_c.keys() if t not in ["cellID", "filename", "n_annotation"]]
    lst_taxa_254 = list(set(curate_biigle_report.DCT_BIIGLE254.values()))
    lst_taxa_839 = [t for t in df_b.keys() if t not in ["cellID", "filename", "image_size_sqpx"] + list(curate_biigle_report.DCT_BIIGLE254.values())]

    df_b_254 = df_b[["cellID", "filename", "image_size_sqpx"] + lst_taxa_254]
    df_b_839 = df_b[["cellID", "filename", "image_size_sqpx"] + lst_taxa_839]

    df_c = df_c.drop(index=df_c[~df_c.filename.isin(df_m[df_m.coralnet == 1].filename)].index)
    df_b_254 = df_b_254.drop(index=df_b_254[~df_b_254.filename.isin(df_m[df_m.biigle254 == 1].filename)].index)
    df_b_839 = df_b_839.drop(index=df_b_839[~df_b_839.filename.isin(df_m[df_m.biigle839 == 1].filename)].index)

    df_c.drop(columns=["filename"], inplace=True)
    df_b_254.drop(columns=["filename"], inplace=True)
    df_b_839.drop(columns=["filename"], inplace=True)

    df_c_pc = df_c.groupby(['cellID'], as_index=False).sum()
    df_b_254_pc = df_b_254.groupby(['cellID'], as_index = False).sum()
    df_b_839_pc = df_b_839.groupby(['cellID'], as_index = False).sum()

    df_c_pc[lst_taxa_coralnet] = df_c_pc[lst_taxa_coralnet].div(df_c_pc['n_annotation'].values, axis=0) * 100
    df_b_254_pc[lst_taxa_254] = df_b_254_pc[lst_taxa_254].div(df_b_254_pc['image_size_sqpx'].values, axis=0) * 100
    df_b_839_pc[lst_taxa_839] = df_b_839_pc[lst_taxa_839].div(df_b_839_pc['image_size_sqpx'].values, axis=0) * 100

    df_c_pc.drop(columns=['n_annotation'], inplace=True)
    df_b_254_pc.drop(columns=['image_size_sqpx'], inplace=True)
    df_b_839_pc.drop(columns=['image_size_sqpx'], inplace=True)

    print("\nRemoving the cells where all CoralNet, BIIGLE254 and BIIGLE839 are present ...")
    lst_cellID_intersection = [c for c in df_c_pc.cellID if c in df_b_254_pc.cellID.to_list() and c in df_b_839_pc.cellID.to_list()]
    df_c_pc.drop(index=df_c_pc[~df_c_pc.cellID.isin(lst_cellID_intersection)].index, inplace=True)
    df_b_254_pc.drop(index=df_b_254_pc[~df_b_254_pc.cellID.isin(lst_cellID_intersection)].index, inplace=True)
    df_b_839_pc.drop(index=df_b_839_pc[~df_b_839_pc.cellID.isin(lst_cellID_intersection)].index, inplace=True)

    print("\nMerging datasets ...")
    df_out = pd.merge(df_c_pc, df_b_254_pc, on="cellID", how="left")
    df_out = pd.merge(df_out, df_b_839_pc, on="cellID", how="left")
    df_out["cellID"] = df_out["cellID"].astype(int)

    print("\nSaving data: {} ...".format(fname_o))
    df_out.to_csv(fname_o, index=False)


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    combine_coralnet_biigle(fname_biigle=args.bfname,
                            fname_coralnet=args.cfname,
                            fname_metadata=args.mfname,
                            folder_area_notfull=args.afolder,
                            fname_sqpx=args.pfname,
                            fname_o=args.ofname)


if __name__ == "__main__":
    main()
