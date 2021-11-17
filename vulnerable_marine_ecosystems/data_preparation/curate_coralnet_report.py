import pyreadr
import argparse


# Example:
#   python data_preparation\curate_coralnet_report.py -i C:\Users\cgros\code\IMAS\ARC_Data\annotation\Circumpolar_Annotation_Data.Rdata -o biodata_step3.csv

DCT_CORALNET = {"BH_BrAnt": "bryozoans_hard_branching_antler-bryozoans",
                "BH_BrHead": "bryozoans_hard_branching_coralhead-bryozoans",
                "BH_BrLeaf": "bryozoans_hard_branching_leaf-bryozoans",
                "BH_Encr": "bryozoans_hard_encrusting-bryozoans",
                "BH_Fenestr": "bryozoans_hard_fenestrate-bryozoans",
                "BH_Lettuce": "bryozoans_hard_lettuce-bryozoans",
                "BH_Mass": "bryozoans_hard_massive-bryozoans",
                "B_Purple": "bryozoans_purple-bryozoans",
                "BS_Dendr": "bryozoans_soft_dendroid-bryozoans",
                "BS_FAN": "bryozoans_soft_fan-bryozoans",
                "BS_Foliac": "bryozoans_soft_foliaceous-bryozoans",
                "HydrC_Br": "hydrocorals_branching-stylasterids",
                "S_Amorph": "sponges_massive_simple-porifera",
                "S_Ball": "sponges_massive_round-porifera",
                "S_Barrel": "sponges_massive_barrel-porifera",
                "S_bead": "sponges_bead-porifera",
                "S_Buried": "sponges_massive_cryptic-porifera",
                "S_Creep": "sponges_crust_creeping_ramose-porifera",
                "S_CupCmplt": "sponges_cup_cup_goblet-porifera",
                "S_CupIncmp": "sponges_cup_incomplete_curled-porifera",
                "S_Disc": "sponges_cup_table_disc-porifera",
                "S_Encr": "sponges_crust_encrusting-porifera",
                "S_Er_Br": "sponges_erect_branching-porifera",
                "S_Er_Lam": "sponges_erect_laminar-porifera",
                "S_Er_Palm": "sponges_erect_palmate-porifera",
                "S_Er_Simp": "sponges_erect_simple-porifera",
                "S_Er_St": "sponges_erect_stalked-porifera",
                "S_Tube": "sponges_hollow_tube_chimney-porifera"
                }


def get_parser():
    parser = argparse.ArgumentParser(add_help=False)

    # MANDATORY ARGUMENTS
    mandatory_args = parser.add_argument_group('MANDATORY ARGUMENTS')
    mandatory_args.add_argument('-i', '--ifname', required=True, type=str,
                                help='RData input filename.')
    mandatory_args.add_argument('-o', '--ofname', required=True, type=str,
                                help='CSV filename output.')

    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                               help='Shows function documentation.')

    return parser


def curate_coralnet_report(fname_i, fname_o):
    # Read data
    df = pyreadr.read_r(fname_i)["cover_images"].reset_index()

    # Get n annotation per images
    df["n_annotation"] = df.sum(axis=1) - df["Unscorable"]

    # Rename column
    df.rename(columns={"rownames": "filename"}, inplace=True)

    # Select columns of interest
    df.drop(columns=[c for c in df.keys() if c not in ["filename", "n_annotation"] + list(DCT_CORALNET.keys())], inplace=True)

    # Rename columns
    df.rename(columns=DCT_CORALNET, inplace=True)
    print(df.head())

    # Saving results
    print("Saving results in: {}...".format(fname_o))
    df.to_csv(fname_o, index=False)


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Run function
    curate_coralnet_report(fname_i=args.ifname,
                           fname_o=args.ofname)


if __name__ == "__main__":
    main()
