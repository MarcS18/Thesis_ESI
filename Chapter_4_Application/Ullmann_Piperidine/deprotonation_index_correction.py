import ast
import os
import csv
import pandas as pd
from tqdm import tqdm
from catsd.parse.Chem import read_mol_file, read_xyz_file

"""
Script to correct the atom indexes for deprotonatable functional groups during molsimplify structure generation. 
"""


def main():
    data_file = '2b_general.csv'
    outfile = 'Deprotonation_Corrected_Indexes.csv'
    # Get list of identifiers from .mol files
    # Get list of all mol files in the structure files folder
    items = os.listdir('./Structure_Files/')
    # Strip .mol to get the identifiers
    identifiers = [str(structure).split('.mol')[0] for structure in items if structure.endswith('.mol')]
    deprotonated_index_list = []
    # Loop over identifiers and extract data
    for ident in tqdm(identifiers):
        # Initiate dict containing deprotonation analysis
        data_dict = {'Ligand': ident,
                     'Old Indexes': 'No Value',
                     'New Indexes': 'No Value',
                     'Deprotonation Type': 'No Value',
                     'Note': 'None'}
        # Get original ligand coordinating atom index as a list
        original_ligand_atoms = get_ligand_atoms(data_file, ident)
        original_ligand_atoms.sort()
        data_dict.update({'Old Indexes': original_ligand_atoms})
        mol_symbs, mol_num, mol_coords = read_mol_file(f'./Structure_Files/{ident}.mol')
        # Get generated CuLX complex from molsimplify
        gen_file_prefix = f'{ident}_CuLX_1'
        gen_xyz_file = f'./Unpacked/{gen_file_prefix}/{gen_file_prefix}/{gen_file_prefix}.xyz'
        try:
            # Parse the xyz file
            xyz_symbs, xyz_num, xyz_coords = read_xyz_file(gen_xyz_file)
            # Get the index of Cu and I in the xyz file
        except:
            continue
        i_xyz_index = xyz_symbs.index('I')
        cu_xyz_index = xyz_symbs.index('Cu')
        # Remove I and Cu from the symbol list
        del xyz_symbs[i_xyz_index]
        del xyz_symbs[cu_xyz_index]
        no_match = []
        no_match_count = 0
        # Find indexes of missing hydrogens
        for index, item in enumerate(xyz_symbs):
            # Check is symbol matches across lists
            if xyz_symbs[index] == mol_symbs[index + no_match_count]:
                # Assign matching symbol to variable
                symbol_match = xyz_symbs[index]
                forward_list = []
                backward_list = []
                orig_forward_list = []
                orig_backward_list = []
                # Propagate along the list
                for symbol in xyz_symbs[index:]:
                    if symbol == symbol_match:
                        forward_list.append(symbol)
                    else:
                        break
                # Get the same lists for the original structure
                for symbol in mol_symbs[(index + no_match_count):]:
                    if symbol == symbol_match:
                        orig_forward_list.append(symbol)
                    else:
                        break
                # print(f'Forward list {forward_list}')
                # print(f'Backward list {backward_list}')
                # print(f'Original Forward list {orig_forward_list}')
                # print(f'Original Backward list {orig_backward_list}')
                # Only a H can be missing
                if xyz_symbs[index] == 'H':
                    if len(forward_list) == len(orig_forward_list):
                        pass
                    else:
                        no_match.append(index)
                        no_match_count += 1
                else:
                    pass
                # Check if lists match across structures
            # If symbol does not match
            elif xyz_symbs[index] != mol_symbs[index + no_match_count]:
                no_match.append(index + no_match_count)
                no_match_count += 1
        # Compare length of no_match list to identify deprotonation
        # Indexes should never be equal, if so call an error
        if not no_match:
            data_dict.update({'Deprotonation Type': 'No Deprotonation'})
            data_dict.update({'New Indexes': original_ligand_atoms})
        # Single Deprotonation if len is one lower
        elif len(no_match) == 1:
            data_dict.update({'Deprotonation Type': 'Single Deprotonation'})
            # Check if missing hydrogen is after the second coordinating atom index
            if no_match[0] > original_ligand_atoms[1]:
                data_dict.update({'New Indexes': original_ligand_atoms})
            # If they don't match check where the missing H is
            else:
                # Hydrogen is after the first coordinating atom
                if original_ligand_atoms[0] < no_match[0] < original_ligand_atoms[1]:
                    new_ligand_atoms = [original_ligand_atoms[0], original_ligand_atoms[1] - 1]
                # Hydrogen is before the first coordinating atom
                elif no_match[0] < original_ligand_atoms[0]:
                    new_ligand_atoms = [original_ligand_atoms[0] - 1, original_ligand_atoms[1] - 1]
                else:
                    new_ligand_atoms = 'Error'
                data_dict.update({'New Indexes': new_ligand_atoms})
                data_dict.update({'Note': 'Possibly affected by molSimplify bug'})
        # Double deprotonation if len is two lower
        elif len(no_match) == 2:
            data_dict.update({'Deprotonation Type': 'Double Deprotonation'})
            # Possibilities:
            #   - Both are before the first coordinating atom
            #   - One is before the first and one between 1 and 2
            #   - One is before the first and one after the second
            #   - Both are between first and second
            #   - One is between first and second, one is after second
            #   - Both are after second
            # Check if both hydrogen indexes are greater than the second coordinating atom index
            if all(x > original_ligand_atoms[1] for x in no_match):
                data_dict.update({'New Indexes': original_ligand_atoms})
            # If they don't match check where the missing Hs are
            else:
                # One is before the first and one between 1 and 2
                if no_match[0] < original_ligand_atoms[0] < no_match[1] < original_ligand_atoms[1]:
                    new_ligand_atoms = [original_ligand_atoms[0] - 1, original_ligand_atoms[1] - 2]
                # One is before the first and one after the second
                elif no_match[0] < original_ligand_atoms[0] and no_match[1] > original_ligand_atoms[1]:
                    new_ligand_atoms = [original_ligand_atoms[0] - 1, original_ligand_atoms[1] - 1]
                # Both are between first and second
                elif original_ligand_atoms[0] < no_match[0] < original_ligand_atoms[1] \
                        and original_ligand_atoms[0] < no_match[1] < original_ligand_atoms[1]:
                    new_ligand_atoms = [original_ligand_atoms[0], original_ligand_atoms[1] - 2]
                # One is between first and second, one is after second
                elif original_ligand_atoms[0] < no_match[0] < original_ligand_atoms[1] < no_match[1]:
                    new_ligand_atoms = [original_ligand_atoms[0], original_ligand_atoms[1] - 1]
                # If both are lower than the first coordinating atom
                elif all(x < original_ligand_atoms[0] for x in no_match):
                    new_ligand_atoms = [original_ligand_atoms[0] - 2, original_ligand_atoms[1] - 2]
                else:
                    new_ligand_atoms = 'Error'
                data_dict.update({'New Indexes': new_ligand_atoms})
                data_dict.update({'Note': 'Possibly affected by molSimplify bug'})
        else:
            print(f'Error on Structure {ident}')
        deprotonated_index_list.append(data_dict)
    write_data(outfile, deprotonated_index_list)
    print('Deprotonation Analysis Complete')


def get_ligand_atoms(file, identifier):
    """
    Get the atoms indexes relating to the ligand coordinating atoms
    :param file: file from crossminer containing the coordinating atom data in column 'Coord Atoms'
    :param identifier: CSD Identifier of the ligand
    :return: List of coordinating atom indexes
    """
    df = pd.read_csv(file, index_col='Identifier')
    return ast.literal_eval(df.loc[identifier, 'Coord Atoms'])


def write_data(file, dict_data):
    """
    Write data to output file
    :param file: output file
    :param dict_data: list of descriptor data
    :return:
    """
    with open(file, 'w+', newline='') as f:
        csv_writer = csv.DictWriter(f, fieldnames=dict_data[0].keys())
        csv_writer.writeheader()
        csv_writer.writerows(dict_data)
    return


if __name__ == '__main__':
    main()
