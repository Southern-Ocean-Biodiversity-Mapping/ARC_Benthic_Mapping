import pyreadr
import argparse


# Example:
#   python data_preparation\03_curate_coralnet_report.py -i ..\..\ARC_Data\annotation\Circumpolar_Annotation_Data_500m.RData -o data\biodata\biodata_step3.csv


DCT_CORALNET = {"Bryozoan_Hard_BrAnt": "bryozoans_hard_branching-bryozoans",
                "Bryozoan_Hard_Branching_Head": "bryozoans_hard_branching-bryozoans",
                "Bryozoan_Hard_Branching_Leaf": "bryozoans_hard_branching-bryozoans",
                "Bryozoan_Hard_Branching_Short": "bryozoans_hard_branching-bryozoans",
                "Bryozoan_Hard_Encr": "bryozoans_hard_encrusting-bryozoans",
                "Bryozoan_Hard_Fenestr": "bryozoans_hard_fenestrate-bryozoans",
                "Bryozoan_Hard_Lettuce": "bryozoans_hard_massive-bryozoans",
                "Bryozoan_Purple": "bryozoans_soft_dendroid-bryozoans",
                "Bryozoan_Soft_Dendroid": "bryozoans_soft_dendroid-bryozoans",
                "Bryozoan_Soft_FAN": "bryozoans_soft_dendroid-bryozoans",
                "Bryozoan_Soft_Foliaceous": "bryozoans_soft_foliaceous-bryozoans",
                "Hydrozoa_Hydrocorals_Branching": "hydrocorals_branching-stylasterids",
                "Sponge_Amorph": "sponges_massive_simple-porifera",
                "Sponge_Amorph_Composite": "sponges_massive_composite-porifera",
                "Sponge_Amorph_Orange": "sponges_massive_composite-porifera",
                "Sponge_Amorph_Large": "sponges_massive_simple-porifera",
                "Sponge_Ball": "sponges_massive_globular-porifera",
                "Sponge_Ball_LongSpicules": "sponges_massive_globular-porifera",
                "Sponge_Barrel_Small": "sponges_cup_tube_sac-porifera",
                "Sponge_Barrel": "sponges_cup_barrel-porifera",
                "Sponge_Buried": "sponges_massive_cryptic-porifera",
                "Sponge_Creep": "sponges_crust_creeping-porifera",
                "Sponge_CupCmplt": "sponges_cup_cup_complete-porifera",
                "Sponge_Disc": "sponges_cup_cup_disc-porifera",
                "Sponge_Encrusting": "sponges_crust_encrusting-porifera",
                "Sponge_Erect_Branching": "sponges_erect_3d-porifera",
                "Sponge_Erect_Laminar": "sponges_erect_2d-porifera",
                "Sponge_Erect_Palmate": "sponges_erect_2d-porifera",
                "Sponge_Erect_Simple": "sponges_erect_1d-porifera",
                "Sponge_Erect_Bottlebrush": "sponges_erect_1d-porifera",
                "Sponge_Erect_Stalked": "sponges_erect_stalked-porifera",
                "Sponge_Tube": "sponges_cup_tube_chimney-porifera",
                "Sponge_Cast": "sponges_erect_1d-porifera",
                #"Sponge_Small": "XX",
                "Bryozoan_Hard_Lettuce_maybe_dead": "bryozoans_hard_massive-bryozoans",
                "Bryozoan_Hard_Branching_Antler_maybe_dead": "bryozoans_hard_branching-bryozoans",
                "Bryozoan_Hard_Branching_Leaf_maybe_dead": "bryozoans_hard_branching-bryozoans"
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

    lst_coralnet_labels = list(df.keys())[1:]
    # Get n annotation per images
    df["n_annotation"] = df[lst_coralnet_labels].sum(axis=1)

    # Rename column
    df.rename(columns={"rownames": "filename"}, inplace=True)

    # Select columns of interest
    df.drop(columns=[c for c in df.keys() if c not in ["filename", "n_annotation"] + list(DCT_CORALNET.keys())], inplace=True)

    # Rename columns
    df.rename(columns=DCT_CORALNET, inplace=True)

    # Sum columns with same column name
    df = df.groupby(lambda x: x, axis=1).sum()

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
