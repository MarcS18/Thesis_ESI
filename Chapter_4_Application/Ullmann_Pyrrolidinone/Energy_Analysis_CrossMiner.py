import catsd.data.globalvars as gbl
from catsd.parse.orca import *
import os
import csv
from tqdm import tqdm

"""
Script to calculate the energy diagram values for the complexes in the mechanism. Calculated from computational
chemistry output files from the transition state workflow. For use with CrossMiner based calculations.
"""


def main():

    # Create list of identifiers from folder names
    # List all folders
    items = [f for f in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), f))]
    # Strip at underscore to get identifier
    identifiers = []
    for folder in items:
        ident = str(folder).split('_')[0]
        identifiers.append(ident)
    # Remove duplicates
    identifiers = list(dict.fromkeys(identifiers))
    # Iterate through identifiers and get correct files
    files = list_dir(os.getcwd())
    # Retain orca output files
    files = [f for f in files if 'orca' in f]
    # Files headers
    headers = ['Ligand', 'Starting materials (CuLI) / kcal/mol', 'Ligand Exchange (CuLpyr) / kcal/mol',
               'Relative Energy (TSOA) / kcal/mol', 'Relative Energy (TSSig) / kcal/mol', 'Products (CuLI) / kcal/mol',
               'Activation Energy (TSOA) / kcal/mol', 'TSSig Activation Energy (TSSig) / kcal/mol']
    outfile = 'Activation_Energies_Pyr.csv'
    # Write headers to output file
    with open(outfile, 'w+') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(headers)
    # For each identifier get all relevant file paths in current directory and subdirectories
    for ident in tqdm(identifiers):
        # Create file list
        comp_files = []
        # Create dictionary of energy data (NaN default for missing data, ligand is ligand identifier)
        energy_data = {'le_e': 'NaN', 'le_corr': 'NaN', 'sm_e': 'NaN', 'sm_corr': 'NaN', 'tsoa_e': 'NaN', 'tsoa_corr':
                       'NaN', 'tssig_e': 'NaN', 'tssig_corr': 'NaN', 'ligand': ident}
        # Add file to list if identifier is present
        for file in files:
            if ident in file:
                comp_files.append(file)
        # Only keep .out files
        comp_files = [f for f in comp_files if f.endswith('.out')]
        comp_files = [f for f in comp_files if not f.endswith('_atom53.out')]
        comp_files = [f for f in comp_files if not f.endswith('.smd.out')]
        # If list is empty no files for identifier present
        if len(comp_files) == 0:
            pass
        # Else extract energy data
        else:
            for file in comp_files:
                # Extract relevant data from .out file based on file name and update dictionary
                if 'CuLpip' in file:
                    if 'energy' in file:
                        energy = get_energy(file)
                        energy_data.update({'le_e': energy})
                    elif 'freq' in file:
                        corr = get_correction(file)
                        energy_data.update({'le_corr': corr})
                    else:
                        pass
                elif 'CuLX' in file:
                    if 'energy' in file:
                        energy = get_energy(file)
                        energy_data.update({'sm_e': energy})
                    elif 'freq' in file:
                        corr = get_correction(file)
                        energy_data.update({'sm_corr': corr})
                    else:
                        pass
                elif 'TSOA' in file:
                    if 'energy' in file:
                        energy = get_energy(file)
                        energy_data.update({'tsoa_e': energy})
                    elif 'freq' in file:
                        corr = get_correction(file)
                        energy_data.update({'tsoa_corr': corr})
                    else:
                        pass
                elif 'TSSig' in file:
                    if 'energy' in file:
                        energy = get_energy(file)
                        energy_data.update({'tssig_e': energy})
                    elif 'freq' in file:
                        corr = get_correction(file)
                        energy_data.update({'tssig_corr': corr})
                    else:
                        pass
                else:
                    pass
            # If dictionary value is NaN skip gibbs free energy calculation and return no value
            # Else Calculate Gibbs Free Energy from Single Point Energy and Correction Factor
            # Then Calculate total free energy of all materials in step
            if energy_data['sm_e'] == 'NaN' or energy_data['sm_e'] is None:
                sm_rel = 'No Value'
                sm = 'No Value'
                prod_rel = 'No Value'
            elif energy_data['sm_corr'] == 'NaN' or energy_data['sm_corr'] is None:
                sm_rel = 'No Value'
                sm = 'No Value'
                prod_rel = 'No Value'
            else:
                sm_g = gibbs_free_energy(energy_data['sm_e'], energy_data['sm_corr'])
                sm = sm_g + gbl.dmf_common['pyrrolidinone']['B97-3c'] + gbl.dmf_common['iodobenzene']['B97-3c'] + \
                    2 * gbl.dmf_common['caesium']['B97-3c'] + gbl.dmf_common['carbonate']['B97-3c']
                product = sm_g + gbl.dmf_common['N_phenyl_pyrrolidinone']['B97-3c'] + \
                    2 * gbl.dmf_common['caesium']['B97-3c'] + gbl.dmf_common['iodide']['B97-3c'] + \
                    gbl.dmf_common['hydrogen_carbonate']['B97-3c']
                sm_rel = relative_gibbs_free_energy(sm, sm)
                prod_rel = relative_gibbs_free_energy(sm, product)
            if energy_data['le_e'] == 'NaN' or energy_data['le_e'] is None:
                le_rel = 'No Value'
            elif energy_data['le_corr'] == 'NaN' or energy_data['le_corr'] is None:
                le_rel = 'No Value'
            else:
                le_g = gibbs_free_energy(energy_data['le_e'], energy_data['le_corr'])
                le = le_g + 2 * gbl.dmf_common['caesium']['B97-3c'] + gbl.dmf_common['iodide']['B97-3c'] + \
                    gbl.dmf_common['iodobenzene']['B97-3c'] + gbl.dmf_common['hydrogen_carbonate']['B97-3c']
                le_rel = relative_gibbs_free_energy(sm, le)
            if energy_data['tsoa_e'] == 'NaN' or energy_data['tsoa_e'] is None:
                tsoa_rel = 'No Value'
            elif energy_data['tsoa_corr'] == 'NaN' or energy_data['tsoa_corr'] is None:
                tsoa_rel = 'No Value'
            else:
                tsoa_g = gibbs_free_energy(energy_data['tsoa_e'], energy_data['tsoa_corr'])
                tsoa = tsoa_g + 2 * gbl.dmf_common['caesium']['B97-3c'] + gbl.dmf_common['iodide']['B97-3c'] + \
                    gbl.dmf_common['hydrogen_carbonate']['B97-3c']
                tsoa_rel = relative_gibbs_free_energy(sm, tsoa)
            if energy_data['tssig_e'] == 'NaN' or energy_data['tssig_e'] is None:
                tssig_rel = 'No Value'
            elif energy_data['tssig_corr'] == 'NaN' or energy_data['tssig_corr'] is None:
                tssig_rel = 'No Value'
            else:
                tssig_g = gibbs_free_energy(energy_data['tssig_e'], energy_data['tssig_corr'])
                tssig = tssig_g + 2 * gbl.dmf_common['caesium']['B97-3c'] + gbl.dmf_common['iodide']['B97-3c'] + \
                    gbl.dmf_common['hydrogen_carbonate']['B97-3c']
                tssig_rel = relative_gibbs_free_energy(sm, tssig)
            # Calculate activation energies from the lowest energy structure
            if sm_rel == 'No Value' or le_rel == 'No Value':
                tsoa_act = 'No Value'
                tssig_act = 'No Value'
            else:
                if float(le_rel) < float(sm_rel):
                    tsoa_act = activation_energy(le_rel, tsoa_rel)
                    tssig_act = activation_energy(le_rel, tssig_rel)
                elif float(le_rel) >= float(sm_rel):
                    tsoa_act = tsoa_rel
                    tssig_act = tssig_rel
                else:
                    tsoa_act = 'No Value'
                    tssig_act = 'No Value'
            # Append relative Gibbs Free Energies to data list
            energies = [energy_data['ligand'], format_value(sm_rel), format_value(le_rel), format_value(tsoa_rel),
                        format_value(tssig_rel), format_value(prod_rel), format_value(tsoa_act), format_value(tssig_act)]
            # Write data to output file
            with open(outfile, 'a') as f:
                csv_writer = csv.writer(f)
                csv_writer.writerow(energies)
    print('Analysis Complete')

    return


def gibbs_free_energy(energy, correction):
    """
    Calculate the Gibbs Free Energy
    :param energy: Single point energy from energy calculation
    :param correction: Correction factor from frequency calculation
    :return: Sum of energy and correction factor
    """
    return energy + correction


def activation_energy(min_struc, ts_struc):
    """
    Calculate the activation energy from the minimum energy structure and transition state energy
    :param min_struc: Gibbs free energy of the minimum structure in the reaction path
    :param ts_struc: Gibbs free energy of the transition state structure
    :return: Activation energy
    """
    if min_struc == 'No Value':
        return 'No Value'
    elif ts_struc == 'No Value':
        return 'No Value'
    else:
        # Calculate difference in energy in kcal/mol (630 conversion factor)
        return ts_struc - min_struc


def relative_gibbs_free_energy(starting_material, comp_struc):
    """
    Calculate the relative Gibbs Free energy in comparison to a zero complex
    :param starting_material: Zero point complex
    :param comp_struc: Complex of interest
    :return: If no value - No Value, If value the difference in energy in kcal/mol
    """
    if comp_struc == 'No Value':
        return 'No Value'
    elif starting_material == 'No Value':
        return 'No Value'
    else:
        # Calculate difference in energy in kcal/mol (630 conversion factor)
        return (comp_struc - starting_material) * 630


def list_dir(directory):
    """
    Get list of all files in the current directory and all sub directories
    :param directory: Directory of interest
    :return: list of files in directory and all subdirectories
    """
    r = []
    for root, dirs, files in os.walk(directory):
        for name in files:
            r.append(os.path.join(root, name))
    return r


def format_value(value):
    """
    Format value based on data type
    :param value: Value to be compared
    :return: string if str, float formatted to 1 decimal place if float
    """
    if type(value) == str:
        return str(value)
    elif type(value) == float:
        return "%.1f" % value


if __name__ == '__main__':
    main()
