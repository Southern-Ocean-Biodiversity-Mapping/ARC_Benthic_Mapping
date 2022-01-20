import os
import zipfile
import argparse
import pandas as pd


# Example:
#   python data_preparation\merge_biigle_reports.py -i C:\Users\cgros\Downloads\292_csv_image_annotation_report(2) -t 839-vme-morpho-taxa.csv,254-catami-mobile-indicator-species.csv -o biodata_step1.csv
# python data_preparation\merge_biigle_reports.py -i C:\Users\cgros\Downloads\292_image_annotation_area_report -o biodata_step1.csv


def get_parser():
    parser = argparse.ArgumentParser(add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-i', '--ifolder', required=True, type=str,
                                help='Input folder containing area BIIGLE reports.')
    mandatory_args.add_argument('-o', '--ofname', required=True, type=str,
                                help='CSV filename output.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def merge_biigle_reports(folder_i, fname_o):
    lst_df = []
    for fname_i in os.listdir(folder_i):
        fname_i = os.path.join(folder_i, fname_i)
        if fname_i.endswith(".xlsx"):
            df = pd.read_excel(fname_i, skiprows=1)
            print("Found {} annotations in {}...".format(len(df), fname_i))
            lst_df.append(df)

    df_stacked = pd.concat(lst_df, axis=0)

    print("\nTotal number of annotations: {}...".format(len(df_stacked)))
    print("Saving result in: {}...".format(fname_o))
    df_stacked.to_csv(fname_o, index=False)


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    merge_biigle_reports(folder_i=args.ifolder,
                         fname_o=args.ofname)


if __name__ == "__main__":
    main()
