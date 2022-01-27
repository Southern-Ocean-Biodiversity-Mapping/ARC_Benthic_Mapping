import os
import argparse
from PIL import Image
import numpy as np
import pandas as pd

# Example:
#   python data_preparation\get_px_size_image.py -i R:\IMAS\Antarctic_Seafloor\Clean_Data_For_Permanent_Storage -o px_size.csv


def get_parser():
    parser = argparse.ArgumentParser(add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-i', '--ifolder', required=True, type=str,
                                help='Folder containing the images.')
    mandatory_args.add_argument('-o', '--ofname', required=True, type=str,
                                help='Output csv filename.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def get_px_size_image(folder_i, fname_o):
    dct_ = {"filename": [], "image_size_sqpx": []}
    for survey_name in os.listdir(folder_i):
        survey_folder = os.path.join(folder_i, survey_name)
        if os.path.isdir(survey_folder) and not survey_name.startswith("Annot") and not survey_name.startswith("z"):
            print(survey_name)
            candidate_folder = [f for f in os.listdir(survey_folder) if f.startswith(survey_name+"_3")]
            if len(candidate_folder) != 1:
                print(os.listdir(survey_folder))
            else:
                image_folder = os.path.join(survey_folder, candidate_folder[0])
                for image_fname in os.listdir(image_folder):
                    image_path = os.path.join(image_folder, image_fname)
                    try:
                        img = Image.open(image_path)
                        img.load()
                        data = np.asarray(img)
                        print(data.shape)
                        width, heigth, _ = data.shape
                        del img, data
                        dct_["filename"].append(image_fname)
                        dct_["image_size_sqpx"].append(width * heigth)
                    except:
                        print("error")
                        print(image_fname)
                        pass

    df = pd.DataFrame.from_dict(dct_)
    df.to_csv(fname_o)


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    get_px_size_image(folder_i=args.ifolder,
                      fname_o=args.ofname)


if __name__ == "__main__":
    main()
