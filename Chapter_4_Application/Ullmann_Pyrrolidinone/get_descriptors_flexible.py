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
    method = "TPSS"
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
    tsoa_outfile = f'Structure_Analysis_Descriptors_TSOA_{method}.csv'
    tssig_outfile = f'Structure_Analysis_Descriptors_TSSIG_{method}.csv'
    data_file = 'Deprotonation_Corrected_Indexes.csv'
    print('Extracting Descriptors')
    tsoa_dict_data = []
    tssig_dict_data = []
    # For each identifier get all relevant file paths in current directory and subdirectories
    for ident in tqdm(identifiers):
        tqdm.write(ident)
        tsoa_descriptors = {'Ligand': ident,
                            'Bite Angle (°)': 'No Value',
                            'Change in Bite Angle (°)': 'No Value',
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
                            'ΔCu-L1 (A)': 'No Value',
                            'ΔCu-L2 (A)': 'No Value',
                            'Cu-I (A)': 'No Value',
                            'Cu-C (A)': 'No Value',
                            'C-I (A)': 'No Value',
                            'Cu-I-C (°)': 'No Value',
                            'I-C-Cu (°)': 'No Value',
                            'C-Cu-I (°)': 'No Value',
                            'Lowdin Charge (Cu)': 'No Value',
                            'Lowdin Charge (C)': 'No Value',
                            'Lowdin Charge (I)': 'No Value',
                            'Lowdin Charge (L1)': 'No Value',
                            'Lowdin Charge (L2)': 'No Value',
                            'Bonded Valence (Cu)': 'No Value',
                            'Bonded Valence (C)': 'No Value',
                            'Bonded Valence (I)': 'No Value',
                            'Bonded Valence (L1)': 'No Value',
                            'Bonded Valence (L2)': 'No Value',
                            'Atomic Population (Cu)': 'No Value',
                            'Atomic Population (C)': 'No Value',
                            'Atomic Population (I)': 'No Value',
                            'Atomic Population (L1)': 'No Value',
                            'Atomic Population (L2)': 'No Value',
                            'Bond Order (Cu-I)': 'No Value',
                            'Bond Order (Cu-C)': 'No Value',
                            'Bond Order (C-I)': 'No Value',
                            'Bond Order (Cu-L1)': 'No Value',
                            'Bond Order (Cu-L2)': 'No Value',
                            'Orbital Charge C(s)': 'No Value',
                            'Orbital Charge C(p)': 'No Value',
                            'Orbital Charge C(pz)': 'No Value',
                            'Orbital Charge C(px)': 'No Value',
                            'Orbital Charge C(py)': 'No Value',
                            'Orbital Charge I(s)': 'No Value',
                            'Orbital Charge I(p)': 'No Value',
                            'Orbital Charge I(pz)': 'No Value',
                            'Orbital Charge I(px)': 'No Value',
                            'Orbital Charge I(py)': 'No Value',
                            'Orbital Charge I(d)': 'No Value',
                            'Orbital Charge I(dxz)': 'No Value',
                            'Orbital Charge I(dyz)': 'No Value',
                            'Orbital Charge I(dxy)': 'No Value',
                            'Orbital Charge I(dz2)': 'No Value',
                            'Orbital Charge I(dx2y2)': 'No Value',
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
                            'Orbital Charge L2(dx2y2)': 'No Value',
                            'Magnitude of the Imaginary Frequency (cm-1)': 'No Value'}
        tssig_descriptors = {'Ligand': ident,
                             'Bite Angle (°)': 'No Value',
                             'Change in Bite Angle (°)': 'No Value',
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
                             'ΔCu-L1 (A)': 'No Value',
                             'ΔCu-L2 (A)': 'No Value',
                             'Cu-I (A)': 'No Value',
                             'Cu-N (A)': 'No Value',
                             'Cu-C (A)': 'No Value',
                             'Cu-O (A)': 'No Value',
                             'C-I (A)': 'No Value',
                             'Amide C-N (A)': 'No Value',
                             'Amide C-O (A)': 'No Value',
                             'N-Cu-I (°)': 'No Value',
                             'Cu-I-C (°)': 'No Value',
                             'I-C-N (°)': 'No Value',
                             'C-N-Cu (°)': 'No Value',
                             'Amide O-C-N (°)': 'No Value',
                             'Lowdin Charge (Cu)': 'No Value',
                             'Lowdin Charge (N)': 'No Value',
                             'Lowdin Charge (C)': 'No Value',
                             'Lowdin Charge (I)': 'No Value',
                             'Lowdin Charge (L1)': 'No Value',
                             'Lowdin Charge (L2)': 'No Value',
                             'Lowdin Charge (Amide C)': 'No Value',
                             'Lowdin Charge (Amide O)': 'No Value',
                             'Bonded Valence (Cu)': 'No Value',
                             'Bonded Valence (N)': 'No Value',
                             'Bonded Valence (C)': 'No Value',
                             'Bonded Valence (I)': 'No Value',
                             'Bonded Valence (L1)': 'No Value',
                             'Bonded Valence (L2)': 'No Value',
                             'Bonded Valence (Amide C)': 'No Value',
                             'Bonded Valence (Amide O)': 'No Value',
                             'Atomic Population (Cu)': 'No Value',
                             'Atomic Population (N)': 'No Value',
                             'Atomic Population (C)': 'No Value',
                             'Atomic Population (I)': 'No Value',
                             'Atomic Population (L1)': 'No Value',
                             'Atomic Population (L2)': 'No Value',
                             'Atomic Population (Amide C)': 'No Value',
                             'Atomic Population (Amide O)': 'No Value',
                             'Bond Order (Cu-I)': 'No Value',
                             'Bond Order (Cu-C)': 'No Value',
                             'Bond Order (C-I)': 'No Value',
                             'Bond Order (Cu-N)': 'No Value',
                             'Bond Order (C-N)': 'No Value',
                             'Bond Order (Cu-L1)': 'No Value',
                             'Bond Order (Cu-L2)': 'No Value',
                             'Bond Order (Cu-O)': 'No Value',
                             'Bond Order (Amide C-O)': 'No Value',
                             'Bond Order (Amide C-N)': 'No Value',
                             'Orbital Charge C(s)': 'No Value',
                             'Orbital Charge C(p)': 'No Value',
                             'Orbital Charge C(pz)': 'No Value',
                             'Orbital Charge C(px)': 'No Value',
                             'Orbital Charge C(py)': 'No Value',
                             'Orbital Charge N(s)': 'No Value',
                             'Orbital Charge N(p)': 'No Value',
                             'Orbital Charge N(pz)': 'No Value',
                             'Orbital Charge N(px)': 'No Value',
                             'Orbital Charge N(py)': 'No Value',
                             'Orbital Charge I(s)': 'No Value',
                             'Orbital Charge I(p)': 'No Value',
                             'Orbital Charge I(pz)': 'No Value',
                             'Orbital Charge I(px)': 'No Value',
                             'Orbital Charge I(py)': 'No Value',
                             'Orbital Charge I(d)': 'No Value',
                             'Orbital Charge I(dxz)': 'No Value',
                             'Orbital Charge I(dyz)': 'No Value',
                             'Orbital Charge I(dxy)': 'No Value',
                             'Orbital Charge I(dz2)': 'No Value',
                             'Orbital Charge I(dx2y2)': 'No Value',
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
                             'Orbital Charge L2(dx2y2)': 'No Value',
                             'Orbital Charge Amide C(s)': 'No Value',
                             'Orbital Charge Amide C(p)': 'No Value',
                             'Orbital Charge Amide C(pz)': 'No Value',
                             'Orbital Charge Amide C(px)': 'No Value',
                             'Orbital Charge Amide C(py)': 'No Value',
                             'Orbital Charge Amide O(s)': 'No Value',
                             'Orbital Charge Amide O(p)': 'No Value',
                             'Orbital Charge Amide O(pz)': 'No Value',
                             'Orbital Charge Amide O(px)': 'No Value',
                             'Orbital Charge Amide O(py)': 'No Value',
                             'Magnitude of the Imaginary Frequency (cm-1)': 'No Value'}
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
        # Calculate original L1-Cu-L2 angle
        original_ligand_angle = calc_angle(x_coordinates, culx_ligand_atoms[0], 0, culx_ligand_atoms[1])
        # Calculate original L1-Cu and Cu-L2 bond distances
        original_l1_distance = calc_distance(x_coordinates, 0, culx_ligand_atoms[0])
        original_l2_distance = calc_distance(x_coordinates, 0, culx_ligand_atoms[1])
        x_freq_file = f'./{x_file_prefix}/{x_file_prefix}_orca_freq.out'
        x_xyz_file = f'./{x_file_prefix}/{x_file_prefix}_orca_opt.xyz'
        # x_atoms, x_num_atoms, x_coords = read_xyz_file(x_xyz_file)
        # x_i_index = x_atoms.index('I')
        # Analyse CuLpip structure
        nu_file_prefix = f'{ident}_CuLpyr_1'
        nu_energy_file = f'./{nu_file_prefix}/{nu_file_prefix}_orca_{method}_energy.out'
        nu_freq_file = f'./{nu_file_prefix}/{nu_file_prefix}_orca_freq.out'
        nu_xyz_file = f'./{nu_file_prefix}/{nu_file_prefix}_orca_opt.xyz'
        # Assuming piperidine was added to the structure last
        culpyr_ligand_atoms = [atom + 1 for atom in coord_atoms]
        # Analyse TSOA structure
        tsoa_file_prefix = f'{ident}_TSOA_pyr_1'
        tsoa_ligand_atoms = [atom + 25 for atom in coord_atoms]
        ts_active_atoms = [0, 4, 13]
        tsoa_bonds = [[4, 13], [0, 4]]
        tsoa_atoms_of_interest = {'c': {'index': 4,
                                        'symbol': 'C',
                                        'orbitals': ['s', 'p', 'pz', 'px', 'py']},
                                  'i': {'index': 13,
                                        'symbol': 'I',
                                        'orbitals': ['s', 'p', 'pz', 'px', 'py', 'd', 'dxz', 'dyz', 'dxy', 'dz2', 'dx2y2']},
                                  'cu': {'index': 0,
                                         'symbol': 'Cu',
                                         'orbitals': ['s', 'p', 'pz', 'px', 'py', 'd', 'dxz', 'dyz', 'dxy', 'dz2', 'dx2y2']},
                                  'l1': {'index': tsoa_ligand_atoms[0],
                                         'symbol': 'L1',
                                         'orbitals': ['s', 'p', 'pz', 'px', 'py', 'd', 'dxz', 'dyz', 'dxy', 'dz2', 'dx2y2']},
                                  'l2': {'index': tsoa_ligand_atoms[1],
                                         'symbol': 'L2',
                                         'orbitals': ['s', 'p', 'pz', 'px', 'py', 'd', 'dxz', 'dyz', 'dxy', 'dz2', 'dx2y2']}
                                  }

        tsoa_energy_file = f'./{tsoa_file_prefix}/{tsoa_file_prefix}_orca_{method}_energy.out'
        tsoa_energy_data = cclib.io.ccread(tsoa_energy_file)
        if tsoa_energy_data.metadata['success']:
            # Get atomic coordinates from file, will fail if the calculation failed
            try:
                tsoa_coordinates = tsoa_energy_data.atomcoords
            except:
                continue
            if len(tsoa_coordinates[0]) == 0:
                continue
            # Check to see if the length of atoms in the ligands is correct (some cases where A H is between two atoms and is deprotonated incorrectly as it can be considered bound to both)
            if len(tsoa_coordinates[0]) == (len(x_coordinates[0]) + 23):
                pass
            else:
                continue
            # Energy of the HOMO
            tsoa_homo_energy = desc.Orca.homo(tsoa_energy_data)
            tsoa_descriptors.update({'HOMO Energy (eV)': '%.4f' % tsoa_homo_energy})
            # Energy of the LUMO
            tsoa_lumo_energy = desc.Orca.lumo(tsoa_energy_data)
            tsoa_descriptors.update({'LUMO Energy (eV)': '%.4f' % tsoa_lumo_energy})
            # Bond Lengths
            tsoa_descriptors.update({'Cu-L1 (A)': '%.2f' % calc_distance(tsoa_coordinates, 0, tsoa_ligand_atoms[0])})
            tsoa_descriptors.update({'Cu-L2 (A)': '%.2f' % calc_distance(tsoa_coordinates, 0, tsoa_ligand_atoms[1])})
            tsoa_descriptors.update({'Cu-I (A)': '%.2f' % calc_distance(tsoa_coordinates, 0, 13)})
            tsoa_descriptors.update({'Cu-C (A)': '%.2f' % calc_distance(tsoa_coordinates, 0, 4)})
            tsoa_descriptors.update({'C-I (A)': '%.2f' % calc_distance(tsoa_coordinates, 4, 13)})
            # Change in bond lengths
            tsoa_descriptors.update({'ΔCu-L1 (A)': '%.2f' % (original_l1_distance - calc_distance(tsoa_coordinates, 0, tsoa_ligand_atoms[0]))})
            tsoa_descriptors.update({'ΔCu-L2 (A)': '%.2f' % (original_l2_distance - calc_distance(tsoa_coordinates, 0, tsoa_ligand_atoms[1]))})
            # Bond Angles
            tsoa_descriptors.update({'Bite Angle (°)': '%.2f' % calc_angle(tsoa_coordinates, tsoa_ligand_atoms[0], 0, tsoa_ligand_atoms[1])})
            tsoa_descriptors.update({'Cu-I-C (°)': '%.2f' % calc_angle(tsoa_coordinates, 0, 13, 4)})
            tsoa_descriptors.update({'I-C-Cu (°)': '%.2f' % calc_angle(tsoa_coordinates, 13, 4, 0)})
            tsoa_descriptors.update({'C-Cu-I (°)': '%.2f' % calc_angle(tsoa_coordinates, 4, 0, 13)})
            # Change in bond angles
            tsoa_descriptors.update({'Change in Bite Angle (°)': '%.2f' % (original_ligand_angle - calc_angle(tsoa_coordinates, tsoa_ligand_atoms[0], 0, tsoa_ligand_atoms[1]))})
            # Lowdin Charges
            tsoa_descriptors.update({'Lowdin Charge (Cu)': '%.4f' % desc.Orca.lowdin_charge(tsoa_energy_data, 0)})
            tsoa_descriptors.update({'Lowdin Charge (C)': '%.4f' % desc.Orca.lowdin_charge(tsoa_energy_data, 4)})
            tsoa_descriptors.update({'Lowdin Charge (I)': '%.4f' % desc.Orca.lowdin_charge(tsoa_energy_data, 13)})
            tsoa_descriptors.update({'Lowdin Charge (L1)': '%.4f' % desc.Orca.lowdin_charge(tsoa_energy_data, tsoa_ligand_atoms[0])})
            tsoa_descriptors.update({'Lowdin Charge (L2)': '%.4f' % desc.Orca.lowdin_charge(tsoa_energy_data, tsoa_ligand_atoms[1])})
            # Bonded Valence
            tsoa_pop_analysis = tsoa_energy_data.mayerpop
            tsoa_descriptors.update({'Bonded Valence (Cu)': '%.4f' % tsoa_pop_analysis[str(0)]['Mayer\'s Bonded Valence']})
            tsoa_descriptors.update({'Bonded Valence (C)': '%.4f' % tsoa_pop_analysis[str(4)]['Mayer\'s Bonded Valence']})
            tsoa_descriptors.update({'Bonded Valence (I)': '%.4f' % tsoa_pop_analysis[str(13)]['Mayer\'s Bonded Valence']})
            tsoa_descriptors.update({'Bonded Valence (L1)': '%.4f' % tsoa_pop_analysis[str(tsoa_ligand_atoms[0])]['Mayer\'s Bonded Valence']})
            tsoa_descriptors.update({'Bonded Valence (L2)': '%.4f' % tsoa_pop_analysis[str(tsoa_ligand_atoms[1])]['Mayer\'s Bonded Valence']})
            # Atomic Population
            tsoa_descriptors.update({'Atomic Population (Cu)': '%.4f' % tsoa_pop_analysis[str(0)]['Mulliken Gross Atomic Population']})
            tsoa_descriptors.update({'Atomic Population (C)': '%.4f' % tsoa_pop_analysis[str(4)]['Mulliken Gross Atomic Population']})
            tsoa_descriptors.update({'Atomic Population (I)': '%.4f' % tsoa_pop_analysis[str(13)]['Mulliken Gross Atomic Population']})
            tsoa_descriptors.update({'Atomic Population (L1)': '%.4f' % tsoa_pop_analysis[str(tsoa_ligand_atoms[0])]['Mulliken Gross Atomic Population']})
            tsoa_descriptors.update({'Atomic Population (L2)': '%.4f' % tsoa_pop_analysis[str(tsoa_ligand_atoms[1])]['Mulliken Gross Atomic Population']})
            # Bond Order
            tsoa_bond_orders = tsoa_energy_data.mayerbo
            tsoa_descriptors.update({'Bond Order (Cu-I)': '%.4f' % get_bond_order(tsoa_bond_orders, 0, 13)})
            tsoa_descriptors.update({'Bond Order (Cu-C)': '%.4f' % get_bond_order(tsoa_bond_orders, 0, 4)})
            tsoa_descriptors.update({'Bond Order (C-I)': '%.4f' % get_bond_order(tsoa_bond_orders, 4, 13)})
            tsoa_descriptors.update({'Bond Order (Cu-L1)': '%.4f' % get_bond_order(tsoa_bond_orders, 0, tsoa_ligand_atoms[0])})
            tsoa_descriptors.update({'Bond Order (Cu-L2)': '%.4f' % get_bond_order(tsoa_bond_orders, 0, tsoa_ligand_atoms[1])})
            # Orbital charges
            tsoa_ocharges = tsoa_energy_data.orbitalcharges
            for atom, values in tsoa_atoms_of_interest.items():
                for orbital in values['orbitals']:
                    get_orbital_charge(tsoa_ocharges, values, orbital, tsoa_descriptors)
            # Midpoint
            tsoa_midpoint = calc_midpoint(tsoa_coordinates, tsoa_ligand_atoms[0], tsoa_ligand_atoms[1])
        # Data from Freq file
        tsoa_freq_file = f'./{tsoa_file_prefix}/{tsoa_file_prefix}_orca_freq.out'
        try:
            tsoa_freq_data = cclib.io.ccread(tsoa_freq_file)
        except:
            continue
        if tsoa_freq_data.metadata['success']:
            # Magnitude of the imaginary frequency
            tsoa_imag_freq = desc.Orca.imaginary_frequency(tsoa_freq_data)
            tsoa_descriptors.update({'Magnitude of the Imaginary Frequency (cm-1)': '%.2f' % tsoa_imag_freq})
            # Get data from the xyz file
            tsoa_xyz_file = f'./{tsoa_file_prefix}/{tsoa_file_prefix}_orca_optts.xyz'
            tsoa_dummy_xyz_file = f'./{tsoa_file_prefix}/{tsoa_file_prefix}_orca_optts_withDummy_NoSubstrates.xyz'
            # Literature Descriptors
            # Remove substrates, add a dummy H at the midpoint between the two designated coordinating atoms and write to a
            # new xyz file
            try:
                tsoa_midpoint
            except NameError:
                tsoa_midpoint = None
            if tsoa_midpoint is not None:
                if os.path.exists(tsoa_dummy_xyz_file):
                    pass
                else:
                    write_dummy_structure(tsoa_xyz_file, tsoa_dummy_xyz_file, tsoa_midpoint)
                # Define new coordinating atom numbers for ligand as module uses 1 as the start point for both the
                # Transition state complex and the complex with the dummy and substrates removed
                complex_tsoa_atoms = [atom + 26 for atom in coord_atoms]
                dummy_tsoa_atoms = [atom + 2 for atom in coord_atoms]
                tsoa_complex_elements, tsoa_complex_coordinates = read_xyz(tsoa_xyz_file)
                tsoa_dummy_elements, tsoa_dummy_coordinates = read_xyz(tsoa_dummy_xyz_file)
                # Buried Volume at 3.5A, 5A and 7A
                try:
                    tsoa_bv3 = BuriedVolume(tsoa_dummy_elements, tsoa_dummy_coordinates, 1, include_hs=True,
                                   excluded_atoms=[len(tsoa_dummy_elements)], radius=3.5, z_axis_atoms=[dummy_tsoa_atoms[0]],
                                   xz_plane_atoms=[dummy_tsoa_atoms[0], dummy_tsoa_atoms[1]])
                    tsoa_descriptors.update({'% Buried Volume (3.5A)': '%.1f' % (tsoa_bv3.fraction_buried_volume * 100)})
                except:
                    pass
                try:
                    tsoa_bv5 = BuriedVolume(tsoa_dummy_elements, tsoa_dummy_coordinates, 1, include_hs=True,
                                   excluded_atoms=[len(tsoa_dummy_elements)], radius=5, z_axis_atoms=[dummy_tsoa_atoms[0]],
                                   xz_plane_atoms=[dummy_tsoa_atoms[0], dummy_tsoa_atoms[1]])
                    tsoa_descriptors.update({'% Buried Volume (5.0A)': '%.1f' % (tsoa_bv5.fraction_buried_volume * 100)})
                except:
                    pass
                try:
                    tsoa_bv7 = BuriedVolume(tsoa_dummy_elements, tsoa_dummy_coordinates, 1, include_hs=True,
                                   excluded_atoms=[len(tsoa_dummy_elements)], radius=7, z_axis_atoms=[dummy_tsoa_atoms[0]],
                                   xz_plane_atoms=[dummy_tsoa_atoms[0], dummy_tsoa_atoms[1]])
                    tsoa_descriptors.update({'% Buried Volume (7.0A)': '%.1f' % (tsoa_bv7.fraction_buried_volume * 100)})
                except:
                    pass
                # It is also possible to plot the steric map
                #tsoa_bv3.plot_steric_map(
                    #filename=f'./{tsoa_file_prefix}/{tsoa_file_prefix}/{tsoa_file_prefix}_orca_optts_stericmap_3_5A.png')
                #tsoa_bv5.plot_steric_map(
                    #filename=f'./{tsoa_file_prefix}/{tsoa_file_prefix}/{tsoa_file_prefix}_orca_optts_stericmap_5A.png')
                #tsoa_bv5.plot_steric_map(
                    #filename=f'./{tsoa_file_prefix}/{tsoa_file_prefix}/{tsoa_file_prefix}_orca_optts_stericmap_7A.png')
                # Cone Angle
                try:
                    tsoa_ca = ConeAngle(tsoa_dummy_elements, tsoa_dummy_coordinates, 1)
                    tsoa_descriptors.update({'Cone Angle (°)': '%.2f' % tsoa_ca.cone_angle})
                except:
                    pass
                # Solvent Accessible Surface Areat
                try:
                    tsoa_sasa = SASA(tsoa_complex_elements, tsoa_complex_coordinates)
                    tsoa_descriptors.update({'SASA (A^2)': '%.2f' % tsoa_sasa.area})
                except:
                    pass
                # Sterimol
                try:
                    tsoa_sterimol = Sterimol(tsoa_complex_elements, tsoa_dummy_coordinates, 1, len(tsoa_dummy_elements))
                    tsoa_descriptors.update({'Sterimol B1': '%.2f' % tsoa_sterimol.B_1_value})
                    tsoa_descriptors.update({'Sterimol B5': '%.2f' % tsoa_sterimol.B_5_value})
                    tsoa_descriptors.update({'Sterimol L': '%.2f' % tsoa_sterimol.L_value})
                except:
                    pass
        tsoa_dict_data.append(tsoa_descriptors)
        # Analyse TSSig Structure
        tssig_file_prefix = f'{ident}_TSSig_pyr_1'
        tssig_energy_file = f'./{tssig_file_prefix}/{tssig_file_prefix}_orca_{method}_energy.out'
        tssig_ligand_atoms = [atom + 25 for atom in coord_atoms]
        tssig_active_atoms = [0, 1, 2, 3]
        tsssig_bond = [[1, 3]]
        tssig_atoms_of_interest = {'c': {'index': 3,
                                         'symbol': 'C',
                                         'orbitals': ['s', 'p', 'pz', 'px', 'py']},
                                   'n': {'index': 1,
                                         'symbol': 'N',
                                         'orbitals': ['s', 'p', 'pz', 'px', 'py']},
                                   'amide c': {'index': 14,
                                         'symbol': 'Amide C',
                                         'orbitals': ['s', 'p', 'pz', 'px', 'py']},
                                   'amide o': {'index': 24,
                                               'symbol': 'Amide O',
                                               'orbitals': ['s', 'p', 'pz', 'px', 'py']},
                                   'i': {'index': 2,
                                         'symbol': 'I',
                                         'orbitals': ['s', 'p', 'pz', 'px', 'py', 'd', 'dxz', 'dyz', 'dxy', 'dz2', 'dx2y2']},
                                   'cu': {'index': 0,
                                          'symbol': 'Cu',
                                          'orbitals': ['s', 'p', 'pz', 'px', 'py', 'd', 'dxz', 'dyz', 'dxy', 'dz2', 'dx2y2']},
                                   'l1': {'index': tssig_ligand_atoms[0],
                                          'symbol': 'L1',
                                          'orbitals': ['s', 'p', 'pz', 'px', 'py', 'd', 'dxz', 'dyz', 'dxy', 'dz2', 'dx2y2']},
                                   'l2': {'index': tssig_ligand_atoms[1],
                                          'symbol': 'L2',
                                          'orbitals': ['s', 'p', 'pz', 'px', 'py', 'd', 'dxz', 'dyz', 'dxy', 'dz2', 'dx2y2']}
                                   }
        tssig_energy_data = cclib.io.ccread(tssig_energy_file)
        if tssig_energy_data.metadata['success']:
            # Get atomic coodinates
            try:
                tssig_coordinates = tssig_energy_data.atomcoords
            except:
                continue
            if len(tssig_coordinates[0]) == 0:
                continue
            # Check to see if the length of atoms in the ligands is correct (some cases where A H is between two atoms and is deprotonated incorrectly as it can be considered bound to both)
            if len(tssig_coordinates[0]) == (len(x_coordinates[0]) + 23):
                pass
            else:
                continue
            # Energy of the HOMO
            tssig_homo_energy = desc.Orca.homo(tssig_energy_data)
            tssig_descriptors.update({'HOMO Energy (eV)': '%.4f' % tssig_homo_energy})
            # Energy of the LUMO
            tssig_lumo_energy = desc.Orca.lumo(tssig_energy_data)
            tssig_descriptors.update({'LUMO Energy (eV)': '%.4f' % tssig_lumo_energy})
            # Bond Lengths
            tssig_descriptors.update({'Cu-L1 (A)': '%.2f' % calc_distance(tssig_coordinates, 0, tssig_ligand_atoms[0])})
            tssig_descriptors.update({'Cu-L2 (A)': '%.2f' % calc_distance(tssig_coordinates, 0, tssig_ligand_atoms[1])})
            tssig_descriptors.update({'Cu-I (A)': '%.2f' % calc_distance(tssig_coordinates, 0, 2)})
            tssig_descriptors.update({'Cu-N (A)': '%.2f' % calc_distance(tssig_coordinates, 0, 1)})
            tssig_descriptors.update({'Cu-C (A)': '%.2f' % calc_distance(tssig_coordinates, 0, 3)})
            tssig_descriptors.update({'C-I (A)': '%.2f' % calc_distance(tssig_coordinates, 3, 2)})
            tssig_descriptors.update({'Cu-O (A)': '%.2f' % calc_distance(tssig_coordinates, 0, 24)})
            tssig_descriptors.update({'Amide C-O (A)': '%.2f' % calc_distance(tssig_coordinates, 14, 24)})
            tssig_descriptors.update({'Amide C-N (A)': '%.2f' % calc_distance(tssig_coordinates, 14, 1)})
            # Change in bond lengths
            tssig_descriptors.update({'ΔCu-L1 (A)': '%.2f' % (original_l1_distance - calc_distance(tssig_coordinates, 0, tssig_ligand_atoms[0]))})
            tssig_descriptors.update({'ΔCu-L2 (A)': '%.2f' % (original_l2_distance - calc_distance(tssig_coordinates, 0, tssig_ligand_atoms[1]))})
            # Bond Angles
            tssig_descriptors.update({'Bite Angle (°)': '%.2f' % calc_angle(tssig_coordinates, tssig_ligand_atoms[0], 0, tssig_ligand_atoms[1])})
            tssig_descriptors.update({'N-Cu-I (°)': '%.2f' % calc_angle(tssig_coordinates, 1, 0, 2)})
            tssig_descriptors.update({'Cu-I-C (°)': '%.2f' % calc_angle(tssig_coordinates, 0, 2, 3)})
            tssig_descriptors.update({'I-C-N (°)': '%.2f' % calc_angle(tssig_coordinates, 2, 3, 1)})
            tssig_descriptors.update({'C-N-Cu (°)': '%.2f' % calc_angle(tssig_coordinates, 3, 1, 0)})
            tssig_descriptors.update({'Amide O-C-N (°)': '%.2f' % calc_angle(tssig_coordinates, 24, 14, 1)})
            # Change in bond angles
            tssig_descriptors.update({'Change in Bite Angle (°)': '%.2f' % (original_ligand_angle - calc_angle(tssig_coordinates, tssig_ligand_atoms[0], 0, tssig_ligand_atoms[1]))})
            # Lowdin Charges
            tssig_descriptors.update({'Lowdin Charge (Cu)': '%.4f' % desc.Orca.lowdin_charge(tssig_energy_data, 0)})
            tssig_descriptors.update({'Lowdin Charge (C)': '%.4f' % desc.Orca.lowdin_charge(tssig_energy_data, 3)})
            tssig_descriptors.update({'Lowdin Charge (I)': '%.4f' % desc.Orca.lowdin_charge(tssig_energy_data, 2)})
            tssig_descriptors.update({'Lowdin Charge (N)': '%.4f' % desc.Orca.lowdin_charge(tssig_energy_data, 1)})
            tssig_descriptors.update({'Lowdin Charge (L1)': '%.4f' % desc.Orca.lowdin_charge(tssig_energy_data, tssig_ligand_atoms[0])})
            tssig_descriptors.update({'Lowdin Charge (L2)': '%.4f' % desc.Orca.lowdin_charge(tssig_energy_data, tssig_ligand_atoms[1])})
            tssig_descriptors.update({'Lowdin Charge (Amide C)': '%.4f' % desc.Orca.lowdin_charge(tssig_energy_data, 14)})
            tssig_descriptors.update({'Lowdin Charge (Amide O)': '%.4f' % desc.Orca.lowdin_charge(tssig_energy_data, 24)})
            # Bonded Valence
            tssig_pop_analysis = tssig_energy_data.mayerpop
            tssig_descriptors.update({'Bonded Valence (Cu)': '%.4f' % tssig_pop_analysis[str(0)]['Mayer\'s Bonded Valence']})
            tssig_descriptors.update({'Bonded Valence (C)': '%.4f' % tssig_pop_analysis[str(3)]['Mayer\'s Bonded Valence']})
            tssig_descriptors.update({'Bonded Valence (I)': '%.4f' % tssig_pop_analysis[str(2)]['Mayer\'s Bonded Valence']})
            tssig_descriptors.update({'Bonded Valence (N)': '%.4f' % tssig_pop_analysis[str(1)]['Mayer\'s Bonded Valence']})
            tssig_descriptors.update({'Bonded Valence (L1)': '%.4f' % tssig_pop_analysis[str(tssig_ligand_atoms[0])]['Mayer\'s Bonded Valence']})
            tssig_descriptors.update({'Bonded Valence (L2)': '%.4f' % tssig_pop_analysis[str(tssig_ligand_atoms[1])]['Mayer\'s Bonded Valence']})
            tssig_descriptors.update({'Bonded Valence (Amide C)': '%.4f' % tssig_pop_analysis[str(14)]['Mayer\'s Bonded Valence']})
            tssig_descriptors.update({'Bonded Valence (Amide O)': '%.4f' % tssig_pop_analysis[str(24)]['Mayer\'s Bonded Valence']})
            # Atomic Population
            tssig_descriptors.update({'Atomic Population (Cu)': '%.4f' % tssig_pop_analysis[str(0)]['Mulliken Gross Atomic Population']})
            tssig_descriptors.update({'Atomic Population (C)': '%.4f' % tssig_pop_analysis[str(3)]['Mulliken Gross Atomic Population']})
            tssig_descriptors.update({'Atomic Population (I)': '%.4f' % tssig_pop_analysis[str(2)]['Mulliken Gross Atomic Population']})
            tssig_descriptors.update({'Atomic Population (N)': '%.4f' % tssig_pop_analysis[str(1)]['Mulliken Gross Atomic Population']})
            tssig_descriptors.update({'Atomic Population (L1)': '%.4f' % tssig_pop_analysis[str(tssig_ligand_atoms[0])]['Mulliken Gross Atomic Population']})
            tssig_descriptors.update({'Atomic Population (L2)': '%.4f' % tssig_pop_analysis[str(tssig_ligand_atoms[1])]['Mulliken Gross Atomic Population']})
            tssig_descriptors.update({'Atomic Population (Amide C)': '%.4f' % tssig_pop_analysis[str(14)]['Mulliken Gross Atomic Population']})
            tssig_descriptors.update({'Atomic Population (Amide O)': '%.4f' % tssig_pop_analysis[str(24)]['Mulliken Gross Atomic Population']})
            # Bond Order
            tssig_bond_orders = tssig_energy_data.mayerbo
            tssig_descriptors.update({'Bond Order (Cu-I)': '%.4f' % get_bond_order(tssig_bond_orders, 0, 2)})
            tssig_descriptors.update({'Bond Order (Cu-N)': '%.4f' % get_bond_order(tssig_bond_orders, 0, 1)})
            tssig_descriptors.update({'Bond Order (Cu-C)': '%.4f' % get_bond_order(tssig_bond_orders, 0, 3)})
            tssig_descriptors.update({'Bond Order (C-I)': '%.4f' % get_bond_order(tssig_bond_orders, 3, 2)})
            tssig_descriptors.update({'Bond Order (Cu-L1)': '%.4f' % get_bond_order(tssig_bond_orders, 0, tssig_ligand_atoms[0])})
            tssig_descriptors.update({'Bond Order (Cu-L2)': '%.4f' % get_bond_order(tssig_bond_orders, 0, tssig_ligand_atoms[1])})
            tssig_descriptors.update({'Bond Order (Amide C-O)': '%.4f' % get_bond_order(tssig_bond_orders, 14, 24)})
            tssig_descriptors.update({'Bond Order (Amide C-N)': '%.4f' % get_bond_order(tssig_bond_orders, 14, 1)})
            # Orbital charges
            tssig_ocharges = tssig_energy_data.orbitalcharges
            for atom, values in tssig_atoms_of_interest.items():
                for orbital in values['orbitals']:
                    get_orbital_charge(tssig_ocharges, values, orbital, tssig_descriptors)
            # Midpoint
            tssig_midpoint = calc_midpoint(tssig_coordinates, tssig_ligand_atoms[0], tssig_ligand_atoms[1])
        # Data from Freq file
        tssig_freq_file = f'./{tssig_file_prefix}/{tssig_file_prefix}_orca_freq.out'
        try:
            tssig_freq_data = cclib.io.ccread(tssig_freq_file)
        except:
            continue
        if tssig_freq_data.metadata['success']:
            # Magnitude of the imaginary frequency
            tssig_imag_freq = desc.Orca.imaginary_frequency(tssig_freq_data)
            tssig_descriptors.update({'Magnitude of the Imaginary Frequency (cm-1)': '%.2f' % tssig_imag_freq})
            # Get data from the xyz file
            tssig_xyz_file = f'./{tssig_file_prefix}/{tssig_file_prefix}_orca_optts.xyz'
            tssig_dummy_xyz_file = f'./{tssig_file_prefix}/{tssig_file_prefix}_orca_optts_withDummy_NoSubstrates.xyz'
            # Literature Descriptors
            # Remove substrates, add a dummy H at the midpoint between the two designated coordinating atoms and write to a
            # new xyz file
            try:
                tssig_midpoint
            except NameError:
                tssig_midpoint = None
            if tssig_midpoint is not None:
                if os.path.isfile(tssig_dummy_xyz_file):
                    pass
                else:
                    write_dummy_structure(tssig_xyz_file, tssig_dummy_xyz_file, tssig_midpoint)
                # Define new coordinating atom numbers for ligand as module uses 1 as the start point for both the
                # Transition state complex and the complex with the dummy and substrates removed
                complex_tssig_atoms = [atom + 26 for atom in coord_atoms]
                dummy_tssig_atoms = [atom + 2 for atom in coord_atoms]
                tssig_complex_elements, tssig_complex_coordinates = read_xyz(tssig_xyz_file)
                tssig_dummy_elements, tssig_dummy_coordinates = read_xyz(tssig_dummy_xyz_file)
                # Buried Volume at 3.5A, 5A and 7A
                try:
                    tssig_bv3 = BuriedVolume(tssig_dummy_elements, tssig_dummy_coordinates, 1, include_hs=True,
                                            excluded_atoms=[len(tssig_dummy_elements)], radius=3.5,
                                            z_axis_atoms=[dummy_tssig_atoms[0]],
                                            xz_plane_atoms=[dummy_tssig_atoms[0], dummy_tssig_atoms[1]])
                    tssig_descriptors.update({'% Buried Volume (3.5A)': '%.1f' % (tssig_bv3.fraction_buried_volume * 100)})
                except:
                    pass
                try:
                    tssig_bv5 = BuriedVolume(tssig_dummy_elements, tssig_dummy_coordinates, 1, include_hs=True,
                                            excluded_atoms=[len(tssig_dummy_elements)], radius=5,
                                            z_axis_atoms=[dummy_tssig_atoms[0]],
                                            xz_plane_atoms=[dummy_tssig_atoms[0], dummy_tssig_atoms[1]])
                    tssig_descriptors.update({'% Buried Volume (5.0A)': '%.1f' % (tssig_bv5.fraction_buried_volume * 100)})
                except:
                    pass
                try:
                    tssig_bv7 = BuriedVolume(tssig_dummy_elements, tssig_dummy_coordinates, 1, include_hs=True,
                                            excluded_atoms=[len(tssig_dummy_elements)], radius=7,
                                            z_axis_atoms=[dummy_tssig_atoms[0]],
                                            xz_plane_atoms=[dummy_tssig_atoms[0], dummy_tssig_atoms[1]])
                    tssig_descriptors.update({'% Buried Volume (7.0A)': '%.1f' % (tssig_bv7.fraction_buried_volume * 100)})
                except:
                    pass
                # It is also possible to plot the steric map
                #tssig_bv3.plot_steric_map(
                    #filename=f'./{tssig_file_prefix}/{tssig_file_prefix}/{tssig_file_prefix}_orca_optts_stericmap_3_5A.png')
                #tssig_bv5.plot_steric_map(
                    #filename=f'./{tssig_file_prefix}/{tssig_file_prefix}/{tssig_file_prefix}_orca_optts_stericmap_5A.png')
                #tssig_bv7.plot_steric_map(
                    #filename=f'./{tssig_file_prefix}/{tssig_file_prefix}/{tssig_file_prefix}_orca_optts_stericmap_7A.png')
                # Cone Angle
                try:
                    tssig_ca = ConeAngle(tssig_dummy_elements, tssig_dummy_coordinates, 1)
                    tssig_descriptors.update({'Cone Angle (°)': '%.2f' % tssig_ca.cone_angle})
                except:
                    pass
                # Solvent Accessible Surface Area
                try:
                    tssig_sasa = SASA(tssig_complex_elements, tssig_complex_coordinates)
                    tssig_descriptors.update({'SASA (A^2)': '%.2f' % tssig_sasa.area})
                except:
                    pass
                # Sterimol
                try:
                    tssig_sterimol = Sterimol(tssig_complex_elements, tssig_dummy_coordinates, 1, len(tssig_dummy_elements))
                    tssig_descriptors.update({'Sterimol B1': '%.2f' % tssig_sterimol.B_1_value})
                    tssig_descriptors.update({'Sterimol B5': '%.2f' % tssig_sterimol.B_5_value})
                    tssig_descriptors.update({'Sterimol L': '%.2f' % tssig_sterimol.L_value})
                except:
                    pass
        tssig_dict_data.append(tssig_descriptors)

    write_data(tsoa_outfile, tsoa_dict_data)
    write_data(tssig_outfile, tssig_dict_data)
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


def write_dummy_structure(file, outfile, midpoint):
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
    new_atoms = atoms[:1] + atoms[25:]
    new_num_atoms = num_atoms - 24
    new_coords = coords[:1] + coords[25:]
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


if __name__ == '__main__':
    main()

