import math
import copy
import argparse
import numpy as np
import pandas as pd


# Example:
#   python data_preparation\curate_biigle_report.py -i biodata_step1.csv -o biodata_step2.csv


def get_parser():
    parser = argparse.ArgumentParser(add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-i', '--ifname', required=True, type=str,
                                help='CSV filename to curate.')
    mandatory_args.add_argument('-o', '--ofname', required=True, type=str,
                                help='CSV filename output.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def curate_biigle_reports(fname_i, fname_o):
    # Read dataframe
    df = pd.read_csv(fname_i)

    # Clean df
    df = df.drop(["annotation_label_id", "label_id", "label_name", "firstname", "lastname", "user_id",
                                  "image_id", "shape_id", "annotation_id"], axis=1)
    df["label_hierarchy"] = df["label_hierarchy"].str.replace(" > ", "_")
    df["width"] = df["attributes"].str.split('"width":').str[1].str.split(',').str[0]
    df["height"] = df["attributes"].str.split('"height":').str[1].str.split(',').str[0].str.split('}').str[0]
    df["area"] = df["attributes"].str.split('"area":').str[1].str.split(',').str[0]

    print("\nComputing annotation area in pixel x pixel...")
    for i_row, row in df.iterrows():
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

        df.loc[i_row, "area_annotation_pix"] = area_annotation_pix
        df.loc[i_row, "area_pix"] = area_image_pix

    # Get image area in m2
    df['area'] = df['area'].astype(float)

    # Clean
    df.rename(columns={"label_hierarchy": "label"}, inplace=True)
    df.drop(["attributes", "shape_name", "points", "width", "height", "image_latitude", "image_longitude"], axis=1, inplace=True)
    print(df.head())

    # Gather annotations of each taxa within each image
    # Init dict
    dct_cover = {"filename": [], "area": [], "area_pix": []}
    for k in df["label"].unique():
        if k not in dct_cover:
            dct_cover[k] = []
    # Fill dict
    for f in df["filename"].unique():
        df_cur = df[df["filename"] == f]
        for k in dct_cover.keys():
            if k not in df["label"].unique():
                dct_cover[k].append(df_cur[k].tolist()[0])
            else:
                df_cur_k = df_cur[df_cur["label"] == k]
                if len(df_cur_k):
                    dct_cover[k].append(df_cur_k["area_annotation_pix"].sum())
                else:
                    dct_cover[k].append(0)
    # Convert to data frame
    df_cover = pd.DataFrame.from_dict(dct_cover)
    print(df_cover.head())

    # Save results
    print("Saving COVER result in: {}...".format(fname_o))
    df_cover.to_csv(fname_o, index=False)


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    curate_biigle_reports(fname_i=args.ifname,
                          fname_o=args.ofname)


if __name__ == "__main__":
    main()
