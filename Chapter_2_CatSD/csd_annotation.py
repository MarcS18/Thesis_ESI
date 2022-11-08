import argparse
import os
from ccdc.io import EntryReader
from ccdc.search import Search


"""
Generates annotation .csv for annotating CatSD.
"""


def parse_arguments():

    parser = argparse.ArgumentParser(
        description='Script to extract annotations for CSD-derived entries')

    parser.add_argument('-o', '--output',
                        help='Output csv file [./csd_annotations.csv]',
                        default='./csd_annotations.csv')
    parser.add_argument('-m', '--maximum_csd_structures',
                        help='Number of entries to include from the current CSD',
                        type=int, default=None)

    return parser.parse_args()


def write_annotations(annotation_csv_path, maximum_csd_structures=None):
    number_csd_structures = 0
    with open(annotation_csv_path, 'w') as writer:
        writer.write('identifier,CSD Refcode, database_name, formula,r factor, is_organic, is_organometallic\n')

        # define search parameters for which entries should be included
        settings = Search.Settings()
        settings.not_polymeric = True
        settings.no_disorder = 'Non-hydrogen'
        settings.has_3d_coordinates = True
        settings.must_have_elements = []
        settings.max_hit_structures = 0
        settings.only_organic = False
        settings.max_r_factor = 10.0
        settings.no_errors = False
        settings.no_powder = False
        settings.no_ions = False
        # settings.must_not_have_elements = ['He', 'Be', 'Ne', 'Al', 'Si', 'Ar', 'Sc', 'Ti', 'V',
        #                                    'Cr', 'Ga', 'Ge', 'As', 'Se', 'Kr', 'Rb', 'Y', 'Zr',
        #                                    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        #                                    'Sn', 'Sb', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr',
        #                                    'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
        #                                    'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir',
        #                                    'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Rn', 'Fr',
        #                                    'Ra', 'Ac', 'Th', 'Pa', 'U']
        settings.no_metals = False
        settings.only_organometallic = False

        for entry in EntryReader('CSD'):

            # filter out entries that don't match the search parameters
            if not settings.test(entry):
                continue

            # filter out entries that don't have a crystal property
            try:
                _ = entry.crystal
            except RuntimeError:
                continue

            # write row to CSV
            refcode = entry.identifier
            formula = entry.crystal.formula
            db_name = entry.database_name
            r_factor = entry.r_factor
            is_org = entry.is_organic
            is_orgmet = entry.is_organometallic
            writer.write('{},{},{},{},{},{},{}\n'.format(
                refcode,
                refcode,
                db_name,
                formula if formula else '',
                r_factor if r_factor else '',
                is_org if is_org else '',
                is_orgmet if is_orgmet else ''))

            # exit if the limit is reached
            number_csd_structures = number_csd_structures + 1
            if maximum_csd_structures is not None\
                    and number_csd_structures >= maximum_csd_structures:
                break


def extract_annotations(annotations_out_path, maximum_csd_structures):
    print('Writing annotations to: ' + annotations_out_path)
    if os.path.exists(annotations_out_path):
        os.remove(annotations_out_path)
    write_annotations(annotations_out_path, maximum_csd_structures)


if __name__ == '__main__':
    args = parse_arguments()
    extract_annotations(args.output, args.maximum_csd_structures)