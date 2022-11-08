# -*- coding: utf-8 -*-
import re
import os
import csv
import cclib
import ast
import catsd.descriptors as desc
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.spatial import distance
from catsd.parse.Chem import read_xyz_file, write_xyz_file
from morfeus import BuriedVolume, ConeAngle, Sterimol, SASA, read_xyz

"""
Generate descriptors from ORCA output files
"""


def main():
    # Create list of identifiers from folder names
    method = 'PBO0'
    print('Listing Directory')
    # List all folders
    items = [f for f in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), f))]
    # Strip at underscore to get identifier
    identifiers = []
    print('Compiling Identifiers')
    for folder in sorted(items):
        ident = str(folder).split('_')[0]
        identifiers.append(ident)
    # Remove duplicates
    identifiers = list(dict.fromkeys(identifiers))
    nu_outfile = f'Structure_Analysis_Descriptors_CuLpyr_{method}.csv'
    data_file = 'Deprotonation_Corrected_Indexes.csv'
    print('Extracting Descriptors')
    nu_dict_data = []
    # For each identifier get all relevant file paths in current directory and subdirectories
    for ident in tqdm(identifiers):
        tqdm.write(ident)
        nu_descriptors = {'Ligand': ident,
                             'Bite Angle (°)': 'No Value',
                             'Cone Angle (°)': 'No Value',
                             'Sterimol B1': 'No Value',
                             'Sterimol B5': 'No Value',
                             'Sterimol L': 'No Value',
                             '% Buried Volume (3.5A)': 'No Value',
                             '% Buried Volume (5.0A)': 'No Value',
                             '% Buried Volume (7.0A)': 'No Value',
                             'SASA (A^2)': 'No Value',
                             'HOMO Energy (eV)': 'No Value',
                             'LUMO Energy (eV)': 'No Value',
                             'Cu-L1 (A)': 'No Value',
                             'Cu-L2 (A)': 'No Value',
                             'Cu-N (A)': 'No Value',
                             'Lowdin Charge (Cu)': 'No Value',
                             'Lowdin Charge (N)': 'No Value',
                             'Lowdin Charge (L1)': 'No Value',
                             'Lowdin Charge (L2)': 'No Value',
                             'Bonded Valence (Cu)': 'No Value',
                             'Bonded Valence (N)': 'No Value',
                             'Bonded Valence (L1)': 'No Value',
                             'Bonded Valence (L2)': 'No Value',
                             'Atomic Population (Cu)': 'No Value',
                             'Atomic Population (N)': 'No Value',
                             'Atomic Population (L1)': 'No Value',
                             'Atomic Population (L2)': 'No Value',
                             'Bond Order (Cu-N)': 'No Value',
                             'Bond Order (Cu-L1)': 'No Value',
                             'Bond Order (Cu-L2)': 'No Value',
                             'Orbital Charge N(s)': 'No Value',
                             'Orbital Charge N(p)': 'No Value',
                             'Orbital Charge N(pz)': 'No Value',
                             'Orbital Charge N(px)': 'No Value',
                             'Orbital Charge N(py)': 'No Value',
                             'Orbital Charge Cu(s)': 'No Value',
                             'Orbital Charge Cu(p)': 'No Value',
                             'Orbital Charge Cu(pz)': 'No Value',
                             'Orbital Charge Cu(px)': 'No Value',
                             'Orbital Charge Cu(py)': 'No Value',
                             'Orbital Charge Cu(d)': 'No Value',
                             'Orbital Charge Cu(dxz)': 'No Value',
                             'Orbital Charge Cu(dyz)': 'No Value',
                             'Orbital Charge Cu(dxy)': 'No Value',
                             'Orbital Charge Cu(dz2)': 'No Value',
                             'Orbital Charge Cu(dx2y2)': 'No Value',
                             'Orbital Charge L1(s)': 'No Value',
                             'Orbital Charge L1(p)': 'No Value',
                             'Orbital Charge L1(pz)': 'No Value',
                             'Orbital Charge L1(px)': 'No Value',
                             'Orbital Charge L1(py)': 'No Value',
                             'Orbital Charge L1(d)': 'No Value',
                             'Orbital Charge L1(dxz)': 'No Value',
                             'Orbital Charge L1(dyz)': 'No Value',
                             'Orbital Charge L1(dxy)': 'No Value',
                             'Orbital Charge L1(dz2)': 'No Value',
                             'Orbital Charge L1(dx2y2)': 'No Value',
                             'Orbital Charge L2(s)': 'No Value',
                             'Orbital Charge L2(p)': 'No Value',
                             'Orbital Charge L2(pz)': 'No Value',
                             'Orbital Charge L2(px)': 'No Value',
                             'Orbital Charge L2(py)': 'No Value',
                             'Orbital Charge L2(d)': 'No Value',
                             'Orbital Charge L2(dxz)': 'No Value',
                             'Orbital Charge L2(dyz)': 'No Value',
                             'Orbital Charge L2(dxy)': 'No Value',
                             'Orbital Charge L2(dz2)': 'No Value',
                             'Orbital Charge L2(dx2y2)': 'No Value'}
        coord_atoms = get_ligand_atoms(data_file, ident)
        # Check for errors in coord_atoms
        if coord_atoms == 'Error' or coord_atoms == 'No Value':
            continue
        else:
            pass
        # Analyse CuLI structure
        x_file_prefix = f'{ident}_CuLX_1'
        x_energy_file = f'./{x_file_prefix}/{x_file_prefix}_orca_{method}_energy.out'
        try:
            x_energy_data = cclib.io.ccread(x_energy_file)
            x_coordinates = x_energy_data.atomcoords
        except:
            continue
        # Check if coordinates are empty
        if len(x_coordinates[0]) == 0:
            continue
        # Get corrected ligand atoms
        culx_ligand_atoms = [atom + 1 for atom in coord_atoms]
        # # Calculate original L1-Cu-L2 angle
        # original_ligand_angle = calc_angle(x_coordinates, culx_ligand_atoms[0], 0, culx_ligand_atoms[1])
        # # Calculate original L1-Cu and Cu-L2 bond distances
        # original_l1_distance = calc_distance(x_coordinates, 0, culx_ligand_atoms[0])
        # original_l2_distance = calc_distance(x_coordinates, 0, culx_ligand_atoms[1])
        x_freq_file = f'./{x_file_prefix}/{x_file_prefix}_orca_freq.out'
        x_xyz_file = f'./{x_file_prefix}/{x_file_prefix}_orca_opt.xyz'
        # x_atoms, x_num_atoms, x_coords = read_xyz_file(x_xyz_file)
        # x_i_index = x_atoms.index('I')
        # Analyse CuLpip structure
        nu_file_prefix = f'{ident}_CuLpyr_1'
        nu_energy_file = f'./{nu_file_prefix}/{nu_file_prefix}_orca_{method}_energy.out'
        nu_xyz_file = f'./{nu_file_prefix}/{nu_file_prefix}_orca_opt.xyz'
        nu_atom_list, nu_num_atoms, nu_coords = read_xyz_file(nu_xyz_file)
        pyr = ['H', 'H', 'H', 'H', 'C', 'H', 'C', 'C', 'O', 'H', 'N', 'C']
        nu_index = check_nu_position(nu_atom_list, pyr)
        if nu_index == 1:
            nu_ligand_atoms = [atom + 13 for atom in coord_atoms]
        else:
            nu_ligand_atoms = [atom + 1 for atom in coord_atoms]
        n_index = nu_index + 10
        culpip_atoms_of_interest = {'n': {'index': n_index,
                                         'symbol': 'N',
                                         'orbitals': ['s', 'p', 'pz', 'px', 'py']},
                                   'cu': {'index': 0,
                                          'symbol': 'Cu',
                                          'orbitals': ['s', 'p', 'pz', 'px', 'py', 'd', 'dxz', 'dyz', 'dxy', 'dz2',
                                                       'dx2y2']},
                                   'l1': {'index': nu_ligand_atoms[0],
                                          'symbol': 'L1',
                                          'orbitals': ['s', 'p', 'pz', 'px', 'py', 'd', 'dxz', 'dyz', 'dxy', 'dz2',
                                                       'dx2y2']},
                                   'l2': {'index': nu_ligand_atoms[1],
                                          'symbol': 'L2',
                                          'orbitals': ['s', 'p', 'pz', 'px', 'py', 'd', 'dxz', 'dyz', 'dxy', 'dz2',
                                                       'dx2y2']}
                                   }
        nu_energy_data = cclib.io.ccread(nu_energy_file)
        if nu_energy_data.metadata['success']:
            # Get atomic coodinates
            try:
                nu_coordinates = nu_energy_data.atomcoords
            except:
                continue
            if len(nu_coordinates[0]) == 0:
                continue
            # Check to see if the length of atoms in the ligands is correct (some cases where A H is between two atoms and is deprotonated incorrectly as it can be considered bound to both)
            if len(nu_coordinates[0]) == (len(x_coordinates[0]) + 11):
                pass
            else:
                continue
            # Energy of the HOMO
            nu_homo_energy = desc.Orca.homo(nu_energy_data)
            nu_descriptors.update({'HOMO Energy (eV)': '%.4f' % nu_homo_energy})
            # Energy of the LUMO
            nu_lumo_energy = desc.Orca.lumo(nu_energy_data)
            nu_descriptors.update({'LUMO Energy (eV)': '%.4f' % nu_lumo_energy})
            # Bond Lengths
            nu_descriptors.update({'Cu-L1 (A)': '%.2f' % calc_distance(nu_coordinates, 0, nu_ligand_atoms[0])})
            nu_descriptors.update({'Cu-L2 (A)': '%.2f' % calc_distance(nu_coordinates, 0, nu_ligand_atoms[1])})
            nu_descriptors.update({'Cu-N (A)': '%.2f' % calc_distance(nu_coordinates, 0, n_index)})
            # Bond Angles
            nu_descriptors.update({'Bite Angle (°)': '%.2f' % calc_angle(nu_coordinates, nu_ligand_atoms[0], 0,
                                                                            nu_ligand_atoms[1])})
            # Lowdin Charges
            nu_descriptors.update({'Lowdin Charge (Cu)': '%.4f' % desc.Orca.lowdin_charge(nu_energy_data, 0)})
            nu_descriptors.update({'Lowdin Charge (N)': '%.4f' % desc.Orca.lowdin_charge(nu_energy_data, n_index)})
            nu_descriptors.update(
                {'Lowdin Charge (L1)': '%.4f' % desc.Orca.lowdin_charge(nu_energy_data, nu_ligand_atoms[0])})
            nu_descriptors.update(
                {'Lowdin Charge (L2)': '%.4f' % desc.Orca.lowdin_charge(nu_energy_data, nu_ligand_atoms[1])})
            # Bonded Valence
            nu_pop_analysis = nu_energy_data.mayerpop
            nu_descriptors.update(
                {'Bonded Valence (Cu)': '%.4f' % nu_pop_analysis[str(0)]['Mayer\'s Bonded Valence']})
            nu_descriptors.update(
                {'Bonded Valence (N)': '%.4f' % nu_pop_analysis[str(n_index)]['Mayer\'s Bonded Valence']})
            nu_descriptors.update({'Bonded Valence (L1)': '%.4f' % nu_pop_analysis[str(nu_ligand_atoms[0])][
                'Mayer\'s Bonded Valence']})
            nu_descriptors.update({'Bonded Valence (L2)': '%.4f' % nu_pop_analysis[str(nu_ligand_atoms[1])][
                'Mayer\'s Bonded Valence']})
            # Atomic Population
            nu_descriptors.update(
                {'Atomic Population (Cu)': '%.4f' % nu_pop_analysis[str(0)]['Mulliken Gross Atomic Population']})
            nu_descriptors.update(
                {'Atomic Population (N)': '%.4f' % nu_pop_analysis[str(n_index)]['Mulliken Gross Atomic Population']})
            nu_descriptors.update({'Atomic Population (L1)': '%.4f' % nu_pop_analysis[str(nu_ligand_atoms[0])][
                'Mulliken Gross Atomic Population']})
            nu_descriptors.update({'Atomic Population (L2)': '%.4f' % nu_pop_analysis[str(nu_ligand_atoms[1])][
                'Mulliken Gross Atomic Population']})
            # Bond Order
            nu_pop_orders = nu_energy_data.mayerbo
            nu_descriptors.update({'Bond Order (Cu-N)': '%.4f' % get_bond_order(nu_pop_orders, 0, n_index)})
            nu_descriptors.update(
                {'Bond Order (Cu-L1)': '%.4f' % get_bond_order(nu_pop_orders, 0, nu_ligand_atoms[0])})
            nu_descriptors.update(
                {'Bond Order (Cu-L2)': '%.4f' % get_bond_order(nu_pop_orders, 0, nu_ligand_atoms[1])})
            # Orbital charges
            nu_ocharges = nu_energy_data.orbitalcharges
            for atom, values in culpip_atoms_of_interest.items():
                for orbital in values['orbitals']:
                    get_orbital_charge(nu_ocharges, values, orbital, nu_descriptors)
            # Midpoint
            nu_midpoint = calc_midpoint(nu_coordinates, nu_ligand_atoms[0], nu_ligand_atoms[1])
        # Data from Freq file
        nu_freq_file = f'./{nu_file_prefix}/{nu_file_prefix}_orca_freq.out'
        try:
            nu_freq_data = cclib.io.ccread(nu_freq_file)
        except:
            continue
        if nu_freq_data.metadata['success']:    
            # Get data from the xyz file
            nu_dummy_xyz_file = f'./{nu_file_prefix}/{nu_file_prefix}_orca_optts_withDummy_NoSubstrates.xyz'
            # Literature Descriptors
            # Remove substrates, add a dummy H at the midpoint between the two designated coordinating atoms and write to a
            # new xyz file
            try:
                nu_midpoint
            except NameError:
                nu_midpoint = None
            if nu_midpoint is not None:
                write_dummy_structure(nu_xyz_file, nu_dummy_xyz_file, nu_midpoint, nu_index)
                # Define new coordinating atom numbers for ligand as module uses 1 as the start point for both the
                # Transition state complex and the complex with the dummy and substrates removed (+1 for Cu +1 for module = +2 total)
                if nu_index == 1:
                    complex_nu_atoms = [atom + 14 for atom in coord_atoms]
                    dummy_nu_atoms = [atom + 2 for atom in coord_atoms]
                else:
                    dummy_nu_atoms = [atom + 2 for atom in coord_atoms]
                    complex_nu_atoms = [atom + 2 for atom in coord_atoms]
                nu_complex_elements, nu_complex_coordinates = read_xyz(nu_xyz_file)
                nu_dummy_elements, nu_dummy_coordinates = read_xyz(nu_dummy_xyz_file)
                # Buried Volume at 3.5A, 5A and 7A
                try:
                    tssig_bv3 = BuriedVolume(nu_dummy_elements, nu_dummy_coordinates, 1, include_hs=True,
                                             excluded_atoms=[len(nu_dummy_elements)], radius=3.5,
                                             z_axis_atoms=[dummy_nu_atoms[0]],
                                             xz_plane_atoms=[dummy_nu_atoms[0], dummy_nu_atoms[1]])
                    nu_descriptors.update(
                        {'% Buried Volume (3.5A)': '%.1f' % (tssig_bv3.fraction_buried_volume * 100)})
                except:
                    pass
                try:
                    tssig_bv5 = BuriedVolume(nu_dummy_elements, nu_dummy_coordinates, 1, include_hs=True,
                                             excluded_atoms=[len(nu_dummy_elements)], radius=5,
                                             z_axis_atoms=[dummy_nu_atoms[0]],
                                             xz_plane_atoms=[dummy_nu_atoms[0], dummy_nu_atoms[1]])
                    nu_descriptors.update(
                        {'% Buried Volume (5.0A)': '%.1f' % (tssig_bv5.fraction_buried_volume * 100)})
                except:
                    pass
                try:
                    tssig_bv7 = BuriedVolume(nu_dummy_elements, nu_dummy_coordinates, 1, include_hs=True,
                                             excluded_atoms=[len(nu_dummy_elements)], radius=7,
                                             z_axis_atoms=[dummy_nu_atoms[0]],
                                             xz_plane_atoms=[dummy_nu_atoms[0], dummy_nu_atoms[1]])
                    nu_descriptors.update(
                        {'% Buried Volume (7.0A)': '%.1f' % (tssig_bv7.fraction_buried_volume * 100)})
                except:
                    pass
                # It is also possible to plot the steric map
                # tssig_bv3.plot_steric_map(
                # filename=f'./{tssig_file_prefix}/{tssig_file_prefix}/{tssig_file_prefix}_orca_optts_stericmap_3_5A.png')
                # tssig_bv5.plot_steric_map(
                # filename=f'./{tssig_file_prefix}/{tssig_file_prefix}/{tssig_file_prefix}_orca_optts_stericmap_5A.png')
                # tssig_bv7.plot_steric_map(
                # filename=f'./{tssig_file_prefix}/{tssig_file_prefix}/{tssig_file_prefix}_orca_optts_stericmap_7A.png')
                # Cone Angle
                try:
                    tssig_ca = ConeAngle(nu_dummy_elements, nu_dummy_coordinates, 1)
                    nu_descriptors.update({'Cone Angle (°)': '%.2f' % tssig_ca.cone_angle})
                except:
                    pass
                # Solvent Accessible Surface Area
                try:
                    tssig_sasa = SASA(nu_complex_elements, nu_complex_coordinates)
                    nu_descriptors.update({'SASA (A^2)': '%.2f' % tssig_sasa.area})
                except:
                    pass
                # Sterimol
                try:
                    tssig_sterimol = Sterimol(nu_complex_elements, nu_dummy_coordinates, 1,
                                              len(nu_dummy_elements))
                    nu_descriptors.update({'Sterimol B1': '%.2f' % tssig_sterimol.B_1_value})
                    nu_descriptors.update({'Sterimol B5': '%.2f' % tssig_sterimol.B_5_value})
                    nu_descriptors.update({'Sterimol L': '%.2f' % tssig_sterimol.L_value})
                except:
                    pass
        nu_dict_data.append(nu_descriptors)

    write_data(nu_outfile, nu_dict_data)
    print('Analysis Complete')


def get_ligand_atoms(file, identifier):
    """
    Get the atoms indexes relating to the ligand coordinating atoms
    :param file: file from crossminer containing the coordinating atom data in column 'Coord Atoms'
    :param identifier: CSD Identifier of the ligand
    :return: List of coordinating atom indexes
    """
    df = pd.read_csv(file, index_col='Ligand')
    try:
        new_indexes = ast.literal_eval(str(df.loc[identifier, 'New Indexes']))
    except:
        new_indexes = 'Error'
    if new_indexes == 'Error' or new_indexes == 'No Value':
        return 'Error'
    else:
        return new_indexes


def get_bond_order(bond_orders, a, b):
    """
    Get bond order for a specific bond from cclib list of mayer's bond orders
    :param bond_orders: list of bond order data as a cclib data.mayerbo list
    :param a: atom 1
    :param b: atom 2
    :return: mayer's bond order
    """
    for bond in bond_orders:
        if (bond[0] == str(a) or bond[1] == str(a)) and (bond[0] == str(b) or bond[1] == str(b)):
            return float(bond[2])
    return 0.00


def get_orbital_charge(data, atom_dict, orbital, descriptor_list, ctype='lowdin'):
    """
    get orbital charge from the dict returned from cclib
    :param data: list dict of orbital charge data
    :param atom_dict: dictionary containing the data relating to each atom to be analysed
    :param orbital: orbital to extract
    :param descriptor_list: which list of descriptors to update
    :param ctype: mulliken or lowdin charge
    :return: None
    """
    try:
        descriptor_list.update({f'Orbital Charge {atom_dict["symbol"]}({orbital})': '%.4f' % float(data[ctype][atom_dict["index"]][orbital])})
    except:
        descriptor_list.update({f'Orbital Charge {atom_dict["symbol"]}({orbital})': 'No Value'})
    return


def calc_distance(coordinates, a, b):
    """
    Calculates the distance between two 3D coordinates
    :param coordinates: XYZ coordinates of the structure
    :param a: atom index of the first atom in the bond
    :param b: atom index of the second atom in the bond
    :return: distance between the two atoms in angstroms
    """
    return distance.euclidean(coordinates[0][a], coordinates[0][b])


def calc_angle(coordinates, a, b, c):
    """
    Calculate the angle between three atoms in space
    :param coordinates: 2D numpy array of coordinates from cclib.io.ccread ORCA coordinates
    :param a: index atom 1
    :param b: index of the centre point of the angle (atom 2)
    :param c: index of atom 3
    :return: Angle in degrees
    """
    ba = coordinates[0][b] - coordinates[0][a]
    bc = coordinates[0][b] - coordinates[0][c]
    angle = np.arccos(np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc)))
    return np.degrees(angle)


def calc_midpoint(coordinates, a, b):
    """
    Calculate the midpoint between two atoms in 3D space
    :param coordinates: coordinates of the structure from cclib
    :param a: index of atom 1
    :param b: index of atom 2
    :return: the midpoint as a list of [x, y, z]
    """
    a_coords = coordinates[0][a]
    b_coords = coordinates[0][b]
    return [(a_coords[0] + b_coords[0]) / 2, (a_coords[1] + b_coords[1]) / 2, (a_coords[2] + b_coords[2]) / 2]


def write_dummy_structure(file, outfile, midpoint, nu_index):
    """
    Writes a new xyz file with the substrates removed and a H dummy atom at the midpoint between the two coordinaing
    atoms
    :param file: input file name
    :param outfile: output file name
    :param midpoint: coordinates of the midpoint between the two coordinating atoms
    :return: None
    """
    # Remove substrates, add dummy and write to xyz (TSOA)
    atoms, num_atoms, coords = read_xyz_file(file)
    if nu_index == 1:
        new_atoms = atoms[:1] + atoms[13:]
        new_coords = coords[:1] + coords[13:]
    else:
        new_atoms = atoms[:-12]
        new_coords = coords[:-12]
    new_num_atoms = len(new_atoms)
    new_coords.append(midpoint)
    new_num_atoms += 1
    new_atoms.append('H')
    write_xyz_file(outfile, new_num_atoms, new_atoms, new_coords)
    return


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


def check_nu_position(nu_atom_list, nucleophile):
    """
    Find location of the nucleophile in the structure
    :param nu_atom_list: list of atomic symbols for the compelex
    :param nucleophile: list of atomic symbols of the nucleophile
    :return: index of the first atom of the nucelphile in the complex
    """
    pip_len = len(nucleophile)
    for i in range(len(nu_atom_list)-pip_len + 1):
        if nucleophile == nu_atom_list[i:i + pip_len]:
            return i
        else:
            pass
    return None


if __name__ == '__main__':
    main()

