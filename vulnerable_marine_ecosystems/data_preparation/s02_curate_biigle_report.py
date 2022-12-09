import argparse
import pandas as pd
import os
import sys
sys.path.insert(1, '../../../biigle_scripts')
from biigle.biigle import Api


# Example:
#   python data_preparation\02_curate_biigle_report.py -i data\biodata\biodata_step1.csv -o data\biodata\biodata_step2.csv

DCT_BIIGLE254 = {"Basketstars": "basket_snake_stars-euryalida",
                 "Basketstar-like": "basket_snake_stars-euryalida",
                 "Crinoid - stalked": "crinoid_stalked-crinoid_stalked",
                 "Urchin - regular pencil": "urchin_regular_pencil-cidaroida",
                 "Polychaete - DF in tube": "polychaete-serpulidae",
                 "Polychaete - Scaleworm": "polychaete-serpulidae",
                 "Polychaete - Scaleworm Purple": "polychaete-serpulidae",
                 "Polychaete - Tubeworm": "polychaete-serpulidae",
                 "Polychaete - other": "polychaete-serpulidae"}
DCT_BIIGLE254_DUPLICATE = {"Brachiopoda": "brachiopods-brachiopoda",
                           "Barnacles": "barnacles-bathylasmatidae"}
EMAIL = "charley.gros@gmail.com"
TOKEN = "SVRTBSUtVcQXZHjkNOxI29Zg2yu0nuhw"


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

    # Remove duplicaetes
    len_w_duplicates = len(df)
    df.drop_duplicates(subset=['label_ids', 'image_id', "annotation_area_sqm", "annotation_width_m"],
                       keep='last', inplace=True)
    print("Removing duplicated annotations: {} ...".format(len_w_duplicates - len(df)))

    # Clean df
    df = df.drop(["annotation_id", "shape_id", "shape_name", "annotation_width_m", "annotation_height_m",
                  "annotation_width_px", "annotation_height_px"], axis=1)

    # Get 839 label tree
    api = Api(EMAIL, TOKEN)
    info_label_tree = api.get('label-trees/839').json()["labels"]
    dct_labels = {}
    lst_all_label = [l["id"] for l in info_label_tree]
    for dct_lab in info_label_tree:
        if dct_lab["parent_id"] in lst_all_label:
            parent_name = [l["name"] for l in info_label_tree if l["id"] == dct_lab["parent_id"]][0]
            dct_labels[int(dct_lab["id"])] = parent_name + "_" + dct_lab["name"]
        else:
            dct_labels[int(dct_lab["id"])] = dct_lab["name"]

    df["label_ids"] = df["label_ids"].apply(lambda x: x if type(x) == int else int(x.split(", ")[-1]))

    df["final_label"] = df.loc[df.label_ids.isin(dct_labels.keys())]["label_ids"].apply(lambda x: dct_labels[x] if x in dct_labels else None)

    print("Combining Jan, Victor and charley annotations ...")
    df.loc[df.label_names.isin(DCT_BIIGLE254.keys()), "final_label"] = df.loc[df.label_names.isin(DCT_BIIGLE254.keys())]["label_names"].apply(lambda x: DCT_BIIGLE254.get(x))
    df["jan_label"] = df.loc[df.label_names.isin(DCT_BIIGLE254_DUPLICATE.keys())]["label_names"].apply(lambda x: DCT_BIIGLE254_DUPLICATE[x] if x in DCT_BIIGLE254_DUPLICATE else None)
    # List of images where Jan annotated a species that I also annotated
    lst_im_jan = df[~df.jan_label.isna()]["image_id"].unique()
    df.reset_index(drop=True, inplace=True)
    v_jan = list(DCT_BIIGLE254_DUPLICATE.values())
    for im_jan in lst_im_jan:
        # Annotations of the current image
        df_im_jan = df[df.image_id == im_jan]
        # Annotations of the current image of the species of interest
        df_im_jan_ = df_im_jan[df_im_jan.final_label.isin(v_jan)]
        # If I have not annotated these species on this image
        if len(df_im_jan_) == 0:
            # Then take into account Jan's annotations
            idx_ = df_im_jan[df_im_jan.jan_label.isin(v_jan)].index
            df.loc[idx_, "final_label"] = df.loc[idx_]["jan_label"]
    df.drop("jan_label", axis=1, inplace=True)

    print("Removing not VME indicator taxa ...")
    df = df[~df.final_label.isna()]
    df.reset_index(drop=True, inplace=True)
    print(df.head())
    print(df.final_label.unique())
    print(" ... resulting with {} annotations ...".format(len(df)))

    # Gather annotations of each taxa within each image
    # Init dict
    dct_cover = {"filename": []}
    for k in df["final_label"].unique():
        dct_cover[k] = []
    # Fill dict
    for f in df["image_filename"].unique():
        df_cur = df[df["image_filename"] == f]
        for k in df["final_label"].unique():
            df_cur_k = df_cur[df_cur["final_label"] == k]
            if len(df_cur_k):
                dct_cover[k].append(df_cur_k["annotation_area_sqpx"].sum())
            else:
                dct_cover[k].append(0)
        dct_cover["filename"].append(f)
    # Convert to data frame
    df_cover = pd.DataFrame.from_dict(dct_cover)

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
