import re
import os
import csv
from tqdm import tqdm
import numpy as np
import pandas as pd
from catsd.parse.orca import *
from catsd.parse import Chem


def main():
    # Create list of identifiers from folder names
    # List all folders
    items = [f for f in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), f))]
    # Strip at underscore to get identifier
    identifiers = []
    print('Compiling Identifiers')
    for folder in items:
        ident = str(folder).split('_')[0]
        identifiers.append(ident)
    # Remove duplicates
    identifiers = list(dict.fromkeys(identifiers))
    # Output file headers
    headers = ['Ligand', 'Starting material (CuLI)', 'Ligand Exchange (CuLpip)',
               'Transition State 1 (TSOA)', 'Transition State 2 (TSSig)']
    outfile = 'Structure_Analysis_Pip.csv'
    # Write headers to output file
    with open(outfile, 'w+') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(headers)
    print('Analysing Structures')
    # For each identifier get all relevant file paths in current directory and subdirectories
    for ident in tqdm(identifiers):
        freq_data = {'pip_min': 'NaN', 'x_min': 'NaN', 'tsoa_ts': 'NaN', 'tssig_ts': 'NaN', 'ligand': ident}
        # Analyse CuLI structure
        struc_path = f'{ident}_CuLX_1'
        x_min = perform_analysis(struc_path, mode='stable', tolerance=-20)
        freq_data.update({'x_min': x_min})
        # Analyse CuLpip structure
        struc_path = f'{ident}_CuLpip_1'
        pip_min = perform_analysis(struc_path, mode='stable', tolerance=-20)
        freq_data.update({'pip_min': pip_min})
        # Analyse TSOA structure
        struc_path = f'{ident}_TSOA_pip_1'
        tsoa_bonds = [[4, 13], [0, 4]]
        freqs, overlap, bond_check = perform_analysis(struc_path, mode='ts', bonds=tsoa_bonds, tolerance=-20, ts_tolerance=-40
                                          , ts_freq_tolerance=0.20)
        # Check if TSSig TS is correct
        if freqs == 'Calc Fail':
            freq_data.update({'tsoa_ts': 'Calc Fail'})
        elif freqs == 'Transition State' and overlap is True and bond_check == 'intermediate':
            freq_data.update({'tsoa_ts': 'Correct TS'})
        elif freqs == 'Not Transition State':
            freq_data.update({'tsoa_ts': 'Incorrect Imaginary Frequencies'})
        else:
            freq_data.update({'tsoa_ts': 'Incorrect TS'})
        # Analyse TSSig Structure
        struc_path = f'{ident}_TSSig_pip_1'
        tsssig_bond = [[1, 3]]
        freqs, overlap, bond_check = perform_analysis(struc_path, mode='ts', bonds=tsssig_bond, tolerance=-20, ts_tolerance=-40
                                          , ts_freq_tolerance=0.33)
        # Check if TSSig TS is correct
        if freqs == 'Calc Fail':
            freq_data.update({'tssig_ts': 'Calc Fail'})
        elif freqs == 'Transition State' and overlap is True and bond_check == 'intermediate':
            freq_data.update({'tssig_ts': 'Correct TS'})
        elif freqs == 'Not Transition State':
            freq_data.update({'tssig_ts': 'Incorrect Imaginary Frequencies'})
        else:
            freq_data.update({'tssig_ts': 'Incorrect TS'})
        # Generate list of data for output file
        struc_data = [freq_data['ligand'], freq_data['x_min'], freq_data['pip_min'], freq_data['tsoa_ts'],
                      freq_data['tssig_ts']]
        # Write data to output file
        with open(outfile, 'a') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow(struc_data)
    print('Analysis Complete')


def perform_analysis(file_prefix, mode='stable', bonds=None, tolerance=0, ts_tolerance=-40, ts_freq_tolerance=0.33):
    """
    perform analysis on the structures
    :param file_prefix: output,file
    :param mode: stable - for stable structures, ts - for transition state
    :param bonds: Bond broken/formed in the transition state
    :param tolerance: tolerance for minimum imaginary frequency
    :param ts_tolerance: tolerance for the minimum size the ts frequency is
    :param ts_freq_tolerance: tolerance for the overlap of the ts_mode and bond stretch
    :return:
    """
    freq_file = f'./{file_prefix}/{file_prefix}/{file_prefix}_orca_freq.out'
    xyz_file = f'./{file_prefix}/{file_prefix}/{file_prefix}_orca_optts.xyz'
    if mode == 'stable':
        try:
            # Check is the structure has no imaginary frequencies
            freq_list = get_frequencies_orca(freq_file)
            if any(float(vib) < tolerance for vib in freq_list[0]):
                return 'Not Minimum'
            else:
                return 'Minimum'
        # If no data return fail
        except:
            return 'Calc Fail'
    elif mode == 'ts':
        try:
            # Retrieve atoms and coordinates from the xyz file
            atoms, no_atoms, xyz_coordinates = Chem.read_xyz_file(xyz_file)
            # Transform into numpy array
            coordinates = np.array([np.array(y) for y in xyz_coordinates])
            # Retrieve frequencies from the freq output file.
            freq_list, vib_matrix = get_frequencies_orca(freq_file, n_atoms=no_atoms, mode='ts')
            if float(freq_list[6]) < ts_tolerance and float(freq_list[7]) > tolerance:
                freq_check = 'Transition State'
            else:
                freq_check = 'Not Transition State'
            # Check the imaginary frequency is along the TS bond mode
            lowest_active_freq = check_imaginary_frequencies(freq_list, vib_matrix, bonds, coordinates
                                                             , no_atoms, overlap_tolerance=ts_freq_tolerance)
            # Check if there is a ts active mode
            if lowest_active_freq != 0:
                active_ts = True
            else:
                active_ts = False
            # Check for intermediate bond length
            bond_check = check_bond_length(bonds[0], coordinates, atoms)
            return freq_check, active_ts, bond_check
        # If no data return fail
        except:
            return 'Calc Fail', 'Calc Fail', 'Calc Fail'
    else:
        pass

    return


def divide_chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


def get_frequencies_orca(out_file, n_atoms=None, mode='min'):
    """
    This function extracts the calculated normal modes from a frequency
    calculation with the corresponding frequencies including the optimized
    geometry
    """
    frequency_list = []
    vibration_matrices = []
    with open(out_file, 'r') as _file:
        # Get section containing frequencies
        data = _file.read()
        vib_data = re.findall(r'VIBRATIONAL FREQUENCIES(.*?)NORMAL MODES', data, re.DOTALL)
        # Join string to a list and split on new lines
        vib_data = "".join(vib_data).split('\n')
        # Remove unneeded lines
        vib_data = vib_data[5:-4]
        for line in vib_data:
            # Remove empty entries from split line on spaces
            split_line = list(filter(None, line.split(' ')))
            # Append wavenumber to vib_freqs
            frequency_list.append(split_line[1])
        if mode == 'ts':
            # Get normal mode data
            mode_data = re.findall(r'NORMAL MODES(.*?)IR SPECTRUM', data, re.DOTALL)
            # Join string to a list and split on new lines
            mode_data = "".join(mode_data).split('\n')
            # Remove unneeded data
            mode_data = mode_data[7:-4]
            mode_factors = []
            for line in mode_data:
                # Remove multiple whitespace and create list of data
                line = " ".join(line.split()).split(" ")
                # Create list of data in each line
                mode_factors.append(line)
            # Create blocks based on number of n-modes
            blocks = list(divide_chunks(mode_factors, (n_atoms*3+1)))
            # Get first 6 active modes
            n_mode = blocks[1]
            # Made columns headers from first entry in list
            modes = n_mode[0]
            # Add entry for Index at the start of the entry
            modes.insert(0, 'Index')
            # Remove headers from data
            mode_atom_data = n_mode[1:]
            # Create dataframe
            df = pd.DataFrame(mode_atom_data, columns=modes)
            # Set index to the Index column
            df = df.set_index('Index')
            # Iterate through each normal_mode (df column)
            for (normal_mode, n_mode_data) in df.iteritems():
                #print(f'Normal Mode:{normal_mode}')
                # Split data frame into atom data [x,y,z]
                atom_data = divide_chunks(n_mode_data, 3)
                list_of_atom_data = []
                # For each atom append [x,y,z] data to list
                for atom in atom_data:
                    eigenvalues = [atom[0], atom[1], atom[2]]
                    list_of_atom_data.append(eigenvalues)
                # Convert list to numpy array
                vib = np.array(list_of_atom_data, dtype='float')
                # Round dataframe to 2 decimal places
                #vib = vib.round(2)
                vibration_matrices.append(vib)
    frequency_list = list(map(float, frequency_list))

    return frequency_list, vibration_matrices


def check_imaginary_frequencies(frequency_list, vibration_matrices,
                                bond_breaking_pairs, coordinates, n_atoms, overlap_tolerance=0.33):
    """
    This function checks imaginary frequencies by projecting them onto each of
    the atom pairs that have bonds being formed or broken.
    """
    bond_matrices = []
    # Get number of imaginary frequencies
    n_imag_freqs = sum(n < 0 for n in frequency_list)
    #print("Number of imaginary frequencies = "+str(n_imag_freqs))
    if n_imag_freqs == 0:
        #print("No imaginary Frequencies")
        lowest_freq_active = None
    else:
        for pair in bond_breaking_pairs:
            # Get atom coordinates
            atom_1_coordinates = coordinates[pair[0], :]
            atom_2_coordinates = coordinates[pair[1], :]
            # Find vector for the transition between atom 1 and atom 2
            transition_direction = atom_2_coordinates - atom_1_coordinates
            # Generate matrix for the transition
            transition_matrix = np.zeros((n_atoms, 3))
            transition_matrix[pair[0], :] = transition_direction
            transition_matrix[pair[1], :] = -transition_direction
            bond_matrices.append(transition_matrix)
        lowest_freq_active = 0
        for i in range(n_imag_freqs):
            if i == 0:
                lowest_freq_active = 0
            #print("transition: "+str(i+1))
            # Get vibration of imaginary mode and convert to a 1D array
            frequency_vector = np.ravel(vibration_matrices[i])
            for count, bond_matrix in enumerate(bond_matrices):
                # Convert bond matrix to 1D array (Mode corresponding to the imaginary frequency)
                transition_vector = np.ravel(bond_matrix)
                # @ matrix multiplication dot product
                overlap = (transition_vector/np.linalg.norm(transition_vector)) @ frequency_vector
                #print(bond_breaking_pairs[count], overlap)
                if abs(overlap) > overlap_tolerance:
                    #print("Vibration along the bond")
                    if i == 0:
                        lowest_freq_active += 1
                else:
                    #print("Vibration not along bond")
                    pass

        #print("Lowest imaginary frequency active along: "+str(lowest_freq_active)+" bonds")
    return lowest_freq_active


def check_bond_length(bond, coordinates, atoms):
    """
    Check if the TS active bond is of an intermediate length
    :param bond: list of atom in the bond
    :param coordinates: atomic coordinates
    :param atoms: symbols for atoms
    :return:
    """
    # Dictionary of the colavent radii for TS active atoms
    covalent_radii = {'C': 0.76, 'N': 0.71, 'I': 1.39, 'Cu': 1.32}
    # Get coordinates and symbols for TS active atoms
    atom1_coords = coordinates[bond[0]]
    atom2_coords = coordinates[bond[1]]
    atom1_symb = atoms[bond[0]]
    atom2_symb = atoms[bond[1]]
    # Bond length can be considered an intermediate if the bond length, rij, between atom i and j is:
    # Get covalent radii of each atom in bond from dictionary
    r_i = covalent_radii[atom1_symb]
    r_j = covalent_radii[atom2_symb]
    # Calculate bond length
    r_ij = np.linalg.norm(atom1_coords-atom2_coords)
    # Calculate bond ratio
    bond_ratio = r_ij / (r_i + r_j)
    if 1.7 >= bond_ratio > 1.0:
        return 'intermediate'
    else:
        return 'not intermediate'


if __name__ == '__main__':
    main()

