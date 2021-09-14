import os
import math
import copy
import argparse
import numpy as np
import pandas as pd


# Example:
#   python curate_biigle_report.py -i ..\..\..\biigle_scripts\20210908_biigle_vme.csv -a C:\Users\cgros\data\image_area -o 20210909_biigle_vme.csv


def get_parser():
    parser = argparse.ArgumentParser(add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-i', '--ifname', required=True, type=str,
                                help='CSV filename to curate.')
    mandatory_args.add_argument('-a', '--afolder', required=True, type=str,
                                help='Folder with area info.')
    mandatory_args.add_argument('-o', '--ofname', required=True, type=str,
                                help='CSV filename output.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def curate_biigle_reports(fname_i, folder_area, fname_o):
    # Read dataframe
    df = pd.read_csv(fname_i)

    # Clean df
    df = df.drop(["annotation_label_id", "label_id", "label_name", "firstname", "lastname", "user_id",
                                  "image_id", "shape_id", "annotation_id"], axis=1)
    df["survey"] = df["filename"].str.split('_').str[0]
    df["label_hierarchy"] = df["label_hierarchy"].str.replace(" > ", "_")
    df["width"] = df["attributes"].str.split('"width":').str[1].str.split(',').str[0]
    df["height"] = df["attributes"].str.split('"height":').str[1].str.split(',').str[0]
    df["area"] = df["attributes"].str.split('"area":').str[1].str.split(',').str[0]
    # Fill area missing values
    df_nan_area_rows = df[df["area"].isnull()]
    lst_nan_area_fname = list(set(df_nan_area_rows["filename"].to_list()))
    lst_nan_area_fname_no_extension = [f.split(".")[0] for f in lst_nan_area_fname]
    lst_nan_area_survey = list(set(df_nan_area_rows["survey"].to_list()))
    for survey in lst_nan_area_survey:
        fname_survey = os.path.join(folder_area, survey+"_area.xlsx")
        # TODO: correct for other surveys. Not available now
        if os.path.isfile(fname_survey) and survey == "PS96":
            df_area = pd.read_excel(fname_survey)
            if survey == "PS96":
                lst_nan_area_fname_no_extension_ps96 = [f.split("__")[1] for f in lst_nan_area_fname_no_extension if f.startswith("PS96")]
                df_area = df_area[df_area["image filename"].isin(lst_nan_area_fname_no_extension_ps96)]
            else:
                df_area = df_area[df_area["image filename"].isin(lst_nan_area_fname)]
            df_area = df_area[["image area in m²", "image filename"]]
            for i_row_area, row_area in df_area.iterrows():
                fname_area, area_area = row_area["image filename"], row_area["image area in m²"]
                for i_row, row in df.iterrows():
                    if fname_area in row["filename"]:
                        df.loc[i_row, "area"] = area_area
    df_nan_area_rows_new = df[df["area"].isnull()]
    print("\nDropping {} rows with missing area info...".format(len(df_nan_area_rows_new)))
    print("\tTODO: ASK JAN...")
    df.drop(df_nan_area_rows_new.index, axis="index", inplace=True)

    df_nan_width_rows = df[df["width"].isnull()]
    print("Dropping {} rows with missing width info...".format(len(df_nan_width_rows)))
    print("\tTODO: FETCH THIS INFO...")
    df.drop(df_nan_width_rows.index, axis="index", inplace=True)

    df_nan_height_rows = df[df["height"].isnull()]
    print("Dropping {} rows with missing height info...".format(len(df_nan_height_rows)))
    print("\tTODO: FETCH THIS INFO...")
    df.drop(df_nan_height_rows.index, axis="index", inplace=True)

    df_nan_longitude_rows = df[df["image_longitude"].isnull()]
    print("Dropping {} rows with missing longitude info...".format(len(df_nan_longitude_rows)))
    print("\tTODO: FETCH THIS INFO...")
    df.drop(df_nan_longitude_rows.index, axis="index", inplace=True)

    df_nan_latitude_rows = df[df["image_latitude"].isnull()]
    print("Dropping {} rows with missing latitude info...".format(len(df_nan_latitude_rows)))
    print("\tTODO: FETCH THIS INFO...")
    df.drop(df_nan_latitude_rows.index, axis="index", inplace=True)

    print("\nComputing annotation area in m2...")
    print(df.head())
    for i_row, row in df.iterrows():
        area_image_m2 = float(row["area"])
        area_image_pix = float(row["width"]) * float(row["height"])
        lst_points = row["points"].split("[")[1].split("]")[0]
        lst_points = lst_points.split(',')
        lst_points = [float(p) for p in lst_points]
        if row["shape_name"] == "Circle":
            area_annotation_pix = np.pi * (lst_points[2] ** 2)

        elif row["shape_name"] == "Rectangle":
            length_one_pix = math.sqrt((lst_points[2] - lst_points[0]) ** 2 + (lst_points[3] - lst_points[1]) ** 2)
            length_two_pix = math.sqrt((lst_points[6] - lst_points[0]) ** 2 + (lst_points[7] - lst_points[1]) ** 2)
            area_annotation_pix = length_one_pix * length_two_pix

        else:
            print("ERROR: Unknown shape_name {} ..." .format(row["shape"]))
            exit()

        area_annotation_m2 = area_annotation_pix * area_image_m2 / area_image_pix
        df.loc[i_row, "area_annotation"] = area_annotation_m2

    df['area'] = df['area'].astype(float)

    print("\n\tTODO: ADD ACQUISITION METHOD...")
    print("\n\tTODO: ADD IMAGE QUALITY...")

    # Clean
    df.rename(columns={"image_longitude": "longitude", "image_latitude": "latitude", "label_hierarchy": "label"}, inplace=True)
    df.drop(["attributes", "shape_name", "points", "width", "height"], axis=1, inplace=True)
    print(df.head())

    dct_ = {"filename": [], "survey": [], "longitude": [], "latitude": [], "area": []}
    for k in df["label"].unique():
        if k not in dct_:
            dct_[k] = []

    dct_count, dct_cover = copy.deepcopy(dct_), copy.deepcopy(dct_)
    for f in df["filename"].unique():
        df_cur = df[df["filename"] == f]
        for k in dct_.keys():
            if k not in df["label"].unique():
                dct_count[k].append(df_cur[k].tolist()[0])
                dct_cover[k].append(df_cur[k].tolist()[0])
            else:
                df_cur_k = df_cur[df_cur["label"] == k]
                if len(df_cur_k):
                    dct_count[k].append(len(df_cur_k))
                    dct_cover[k].append(df_cur_k["area_annotation"].sum())
                else:
                    dct_count[k].append(0)
                    dct_cover[k].append(0)

    df_count = pd.DataFrame.from_dict(dct_count)
    df_cover = pd.DataFrame.from_dict(dct_cover)
    print(df_count.head())
    print(df_cover.head())

    print("\n\tTODO: GET EMPTY IMAGES...")

    fname_o_count = fname_o.split(".csv")[0] + "_count.csv"
    print("Saving COUNT result in: {}...".format(fname_o_count))
    df_count.to_csv(fname_o_count, index=False)
    fname_o_cover = fname_o.split(".csv")[0] + "_cover.csv"
    print("Saving COVER result in: {}...".format(fname_o_cover))
    df_count.to_csv(fname_o_cover, index=False)


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    curate_biigle_reports(fname_i=args.ifname, folder_area=args.afolder, fname_o=args.ofname)


if __name__ == "__main__":
    main()
