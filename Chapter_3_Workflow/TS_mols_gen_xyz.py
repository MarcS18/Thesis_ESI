import csv
import os
import pathlib
import shutil
import glob
import argparse
import catsd.data.globalvars as gbl
import catsd.parse.molsimplify as ms

"""
Generates molsimplify input files. For use with 3D .xyz/.mol structure files.
"""


def main():
    # Get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("i", help="Input file name", type=str)
    parser.add_argument("-n", "--nucleophile", help="Key for nucleophile", type=str, choices=["pip", "pyr"], default="pip")
    parser.add_argument("-x", "--halide", help="Key for halide", type=str, choices=["iodide", "bromide", "chloride"],
                        default="iodide")
    args = parser.parse_args()
    infile = args.i
    if args.nucleophile:
        nu = args.nucleophile
    else:
        nu = 'pip'
    # Define variables from command line inputs
    if args.halide:
        x = args.halide
    else:
        x = 'iodide'
    if nu == 'pip':
        ts_list = ['TSOA_pip', 'TSSig_pip']
    elif nu == 'pyr':
        ts_list = ['TSOA_pyr', 'TSSig_pyr']
    else:
        print('Nucleophile not supported')

    csv_data = csv_reader(infile)  # Read input file
    for ligand in csv_data:
        csv_data_dent = input_gen(ligand, nu)  # Generate data from input file
        csv_data_dent_x = input_gen(ligand, x, halide_flag=True)  # Generate halide data
        csv_data_dent_ts = ts_input_gen(ligand)  # Generate transition state data
        # Write molSimplify input files
        input_writer(csv_data_dent)
        input_writer(csv_data_dent_x)
        # Write molSimplify input files for each transition state
        for ts in ts_list:
            ts_input_writer(csv_data_dent_ts, ts)
    # Make output directory and move generated files into it
    out_dir = 'molsimplify_input_files'
    if not os.path.exists(out_dir):
        os.makedirs('molsimplify_input_files')
    input_files = glob.glob(str('*.inp'))
    for file in input_files:
        shutil.move(file, out_dir)

    print('Structure creation complete')

    return()


def input_writer(ligand_data):

    # Write molsimplify input file from catsd/parse/molsimplify instances
    for ligand in ligand_data:
        with open(str(ligand[1]+'.inp'), 'w+') as f:
            ms.complex_name(str(ligand[1]), f)
            ms.core('copper', f)
            ms.oxstate('I', f)
            ms.coord(ligand[4], f)
            ms.geometry(ligand[5], f)
            ms.lig(ligand[0]+', '+ligand[6], f)
            ms.lig_freq(ligand[8], f)
            ms.spin('1', f)
            ms.ff('uff', f)
            ms.ffoption('ba', f)
            ms.keepHs(ligand[10], f)
            ms.ligalign('true', f)
            ms.skipML('true', f)

    return


def ts_input_writer(ts_data, ts_core):

    # Write molsimplify input file for each transition state
    for ligand in ts_data:
        with open(str(ligand[0]+'_'+str(ts_core)+'_'+ligand[1]+'.inp'), 'w+') as f:
            ms.complex_name(str(ligand[0]+'_'+str(ts_core)+'_'+ligand[1]), f)
            ms.core(gbl.ts_cores[ts_core]['core_name'], f)
            ms.oxstate('0', f)  # charge of the core
            ms.spin('1', f)  # spin of the core
            ms.replig('true', f)
            ms.lig(ligand[0], f)
            ms.lig_freq(ligand[7], f)
            ms.ccatoms(gbl.ts_cores[ts_core]['ccatoms'], f)
            ms.ligloc('true', f)
            ms.ligalign('true', f)
            ms.keepHs('auto', f)
            ms.ffoption('c', f)
            ms.skipML('true', f)

    return


def csv_reader(infile):

    # Read csv file and extract number, name, smiles, coordinating atoms and keepHs if applicable

    with open(infile, 'r') as csv_file:
        csv_data = []
        ligand_reader = csv.reader(csv_file)  # Open csv.reader instance
        next(ligand_reader, None)  # Skip header line
        for row in ligand_reader:
            l_id = row[0]  # Ligand ID
            l_name = row[0]
            l_mol = row[0]
            l_index = row[1]  # Ligand index (coordination mode)
            c_atoms = row[4]  # Ligand coordinating atoms
            lig_data = [l_id, l_index, l_name, l_mol, c_atoms]  # Create list from data
            csv_data.append(lig_data)  # Append list of data to csv list to generate entire list from the csv

    return csv_data


def ts_input_gen(ligand):
    """
    Generate list of data for the generation of transition state structures via molsimplify
    :param ligand: ligand data from the csv file
    :return:  list of transition state data
    """
    list_of_ligand_data = []
    l_id = ligand[0] # Ligand ID
    l_coord_atoms = ligand[4]  # Coordinating atoms in the ligand
    l_index = ligand[1]  # Ligand index (coordination mode of ligand)
    l_name = ligand[2]  # Ligand name
    l_mol = ligand[3]  # Ligand mol file

    # Create ligand data for all transition states
    coord_num = '5'  # Coordination number for the generated complex
    coord_geom = 'tbp'  # Coordination geometry for the generated complex
    ligfreq = '1'  # Frequency of the ligand (Default 1 for bidentate ligands)
    ligand_data = [l_id, l_index, l_name, coord_num, coord_geom, str(l_mol), str(l_name), ligfreq, str(l_coord_atoms)]
    list_of_ligand_data.append(ligand_data)

    return list_of_ligand_data


def input_gen(ligand, nu, nu_keeph='False', lig_keeph='auto', halide_flag=False):
    """
    Generate list of data for the generation of non-TS structures via molsimplify
    :param ligand: list of ligand data from the csv file
    :param nu: nucleophile name
    :param nu_keeph: whether to keep hydrogens on the nucleophile (true, false, auto)
    :param lig_keeph: whether to keep hydrogens on the ligand (true, false auto)
    :param halide_flag: whether to use a halide instead of a nucleophile
    :return: list of data for generation of the input file
    """

    list_of_ligand_data = []

    l_id = ligand[0]  # Ligand ID
    l_index = ligand[1]  # Ligand index
    l_name = ligand[2]  # Ligand name

    # Create ligand data for trigonal planar complexes.
    coord_num = '3'  # Coordination number of the complex
    coord_geom = 'tpl'  # Coordination geometry of the complex
    # Data is taken from a dictionary contained in catsd/data/globalvars.py
    # If you are building with a halide
    if halide_flag:
        threed_input = str(gbl.halides[nu]['name'])  # 3D structure input (name in molsimplify database)
        ligatoms = str(gbl.halides[nu]['catom'])  # ligand coordinating atoms (from molsimplify database)
        lignames = str(l_name + ', ' + gbl.halides[nu]['name'])  # Names of the ligands
        complex_name = str(l_id+'_CuLX_'+l_index)  # Name of the output complex
    # If you are building with a supported nucleophile
    else:
        threed_input = str(gbl.nucleophiles[nu]['name'])  # 3D structure input (name in molsimplify database)
        ligatoms = str(gbl.nucleophiles[nu]['catom'])  # ligand coordinating atoms (from molsimplify database)
        lignames = str(l_name + ', ' + gbl.nucleophiles[nu]['name'])  # Names of the ligands
        complex_name = str(l_id+'_CuL'+nu+'_'+l_index)  # Name of the output complex
    ligkeephs = str(lig_keeph + ', ' + nu_keeph)  # Write keepH line from keepHs inputs
    ligfreq = '1, 1'  # Write freq line for the frequency of each ligand
    ligand_data = [l_id, complex_name, l_index, l_name, coord_num, coord_geom, threed_input, lignames, ligfreq, ligatoms, ligkeephs]
    list_of_ligand_data.append(ligand_data)

    return list_of_ligand_data


if __name__ == "__main__":
    main()
