from catsd.data import globalvars as gbl
from catsd.exceptions import XYZfileWrongFormat, XYZfileDidNotExist
import os


def FromXYZ(file):

    if not os.path.exists(file):
        raise XYZfileDidNotExist(f'{file} did not exist')

    if not file.endswith('.xyz'):
        raise XYZfileWrongFormat

    with open(file, 'r') as xyz_file:
        try:
            # First item in an xyz file is the number of atoms
            n_atoms = int(xyz_file.readline().split()[0])

        except (IndexError, ValueError):
            raise XYZfileWrongFormat

        xyz_coords = []
        for line in xyz_file:
            xyz_coords.append(line)
        list_of_coords = xyz_coords[2:]

    return list_of_coords


def ToXYZ(coords, outname, titleline='', append=False):
    """
    Take xyz coordinates as a list and writes an xyz file
    :param coords: List of xyz coordinates
    :param outname: Output file names
    :param titleline: title line of the xyz file
    :param append: Append molecular coordinates to existing file
    :return:
    """
    assert coords is not None
    with open(outname + '.xyz', 'a' if append else 'w') as xyz_file:
        print(len(coords), titleline, sep='\n', file=xyz_file)

        for atom in range(len(coords)):
            print(coords[atom], sep='\n', file=xyz_file)

    return None


def find_metal_index(coords):
    """
    Find metal atoms and return the index in the xyz file as a list
    :param coords: xyz coordinates as a list
    :return: list of metal atoms as xyz index
    """
    for atom in coords:
        atom_idxs = []
        if atom[0] in gbl.metalslist:
            metal_idx = coords.index(atom)
            atom_idxs.append(metal_idx)
            return atom_idxs


def spinmult_to_elec(spin_mult):
    """
    Convert spin multiplicity to number of unpaired electrons
    :param spin_mult: Spin multiplicity
    :return: Number of unpaired electrons
    """
    unpaired_electrons = int(spin_mult) - 1
    return unpaired_electrons


def spinmult_to_spin(spin_mult):
    """
    Convert spin multiplicity to absolute spin
    :param spin_mult: Spin multiplicity
    :return: Absolute spin
    """
    absolute_spin = (int(spin_mult) - 1) / 2
    return absolute_spin


def read_xyz_file(filename):
    """
    read xyz file
    :param filename: Name of the .xyz file
    :return: list of atomic symbols, number of atoms in the structure, list of coordinates
    """

    atomic_symbols = []
    xyz_coordinates = []
    title = ""

    with open(filename, "r") as file:
        for line_number, line in enumerate(file):
            if line_number == 0:
                num_atoms = int(line)
            elif line_number == 1:
                title = line
                if "charge=" in line:
                    charge = int(line.split("=")[1])
            else:
                atomic_symbol, x, y, z = line.split()
                atomic_symbols.append(atomic_symbol)
                xyz_coordinates.append([float(x), float(y), float(z)])

    atoms = [atom for atom in atomic_symbols]

    return atoms, num_atoms, xyz_coordinates


def write_xyz_file(filename, num_atoms, symbols, xyz_coordinates):
    """
    Write an xyz file, to be used with data from read_xyz_file
    :param filename: filename
    :param num_atoms: number of atoms in the structure
    :param symbols: list of atomic symbols for the atoms present
    :param xyz_coordinates: list of xyz coordinates of each atom
    :return:
    """
    with open(filename, "w+") as file:
        file.write(f'{num_atoms}\n\n')
        for i in range(len(symbols)):
            file.write(f'{symbols[i]} {xyz_coordinates[i][0]} {xyz_coordinates[i][1]} {xyz_coordinates[i][2]}\n')
        file.write('\n')
    return


def read_mol_file(filename):
    """
    read mol file
    :param filename: Name of the .mol file
    :return: list of atomic symbols, number of atoms in the structure, list of coordinates
    """
    atomic_symbols = []
    mol_coordinates = []
    bond_data = []
    with open(filename, "r") as file:
        for line_number, line in enumerate(file):
            if line_number == 0:
                title = line
            elif line_number == 1:
                software = line
            elif line == 2:
                comment = line
            elif line_number == 3:
                atom_num = int(float(line[0:3]))
            elif 4 <= line_number <= (3 + atom_num):
                x, y, z, atomic_symbol, a = line.split()
                atomic_symbols.append(atomic_symbol)
                mol_coordinates.append([float(x), float(y), float(z)])
            else:
                pass

    atoms = [atom for atom in atomic_symbols]

    return atoms, atom_num, mol_coordinates
