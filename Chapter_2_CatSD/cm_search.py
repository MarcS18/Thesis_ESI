import argparse
import os
import time
import csv
import glob
import shutil
import math
import re
from decimal import *
from ccdc.pharmacophore import Pharmacophore
from ccdc.io import MoleculeWriter
from ccdc.io import EntryReader
from tqdm import tqdm

"""
Performs a CrossMiner search using the CSD python API. For use with CatSD. Includes identification of coordinating atoms,
duplicate removal, filtering by elements and molecular weight and cleaning up structures.
Saves all structures as a .mol file. Generates ligand .dict for use with structure generation with molSimplify.
Generates input file for TS_mols_gen_xyz.py.
"""

# Set rounding rule to always round down for decimal numbers
getcontext().rounding = ROUND_DOWN


def write_csv(outfile, data):
    """
    Write a csv file containing csd identifier, index, name, xyz/mol file name, coordinating atoms, frequency and rsmd
    :param outfile: output file name
    :param data: data returned from run_search
    :return: None
    """
    with open(outfile+'.csv', 'w+') as of:
        csv_writer = csv.writer(of)
        csv_writer.writerow(['Identifier', 'Index', 'Chemical Name', 'Structure File', 'Coord Atoms', 'Freq', 'rsmd'])
        for line in data:
            csv_writer.writerow(line)
    return


def get_point_coords(point):
    """
    Get xyz coordinates from a feature point
    :param point: str(feature_point)
    :return: xyz coordinates
    """
    feat_point = str(point)
    # Strip surrounding brackets
    feat_point = feat_point[1:-1]
    # Get numbers between brackets
    coords = feat_point[feat_point.find('(')+1:feat_point.find(')')]
    # Split coords, format and rejoin
    point_coords = coords.split(',')
    coords = ', '.join('%.1f' % round(Decimal(float(c)), 1) for c in point_coords)
    coords = [float(x) for x in coords.split(',')]
    # Create inverted set with flipped sign to account for inverted coordinates
    inverted_coords = ', '.join('%.1f' % round(Decimal(float(c)) * -1, 1) for c in point_coords)
    inverted_coords = [float(x) for x in inverted_coords.split(',')]
    set_of_coords = [coords, inverted_coords]
    return set_of_coords


class Runner(argparse.ArgumentParser):
    """
    Runs the searches
    """

    def __init__(self):
        """
        Arguments for running the script.
        """
        super(self.__class__, self).__init__(description=__doc__)
        self.add_argument(
            '-n', '--name', default='unknown',
            help='Name of the search'
        )
        self.add_argument(
            '-d', '--database', default=None,
            help='Feature database to search'
        )
        self.add_argument(
            '-m', '--max-hit-structures', type=int, default=50000,
            help='Maximum number of results to return.'
        )
        self.add_argument(
            '-c', '--catalophore', default=None,
            help='Catalophore to use for the search'
        )
        self.add_argument(
            '-t', '--threads', type=int, default=4,
            help='Number of threads'
        )
        self.add_argument(
            '-r', '--rmsd', default=1,
            help='RMSD used for matching structures'
        )
        self.add_argument(
            '-w', '--max-molecular-weight', type=int, default=500,
            help='Maximum molecular weight of hit structures'
        )
        self.args = self.parse_args()

    def run(self):
        """
        Run script
        :return:
        """
        out_dir = 'Structure_Files_' + self.args.name
        # Load database. default csd database if no input, else load user-specified feature database.
        print('Loading Feature Database')
        if self.args.database is None:
            self.database = Pharmacophore.default_feature_database()
        else:
            self.database = Pharmacophore.FeatureDatabase.from_file(self.args.database)
        # Run search
        list_of_data = self.run_search()
        # Write csv file for generation of molsimplify input files
        write_csv(self.args.name, list_of_data)
        # Make output directory for structure files
        if not os.path.exists(out_dir):
            os.makedirs('Structure_Files_' + self.args.name)
        # Glob all mol files from search
        input_files = glob.glob(str('*.mol'))
        # Move mol files into structure file directory
        for file in input_files:
            shutil.move(file, out_dir)
        # Move ligand.dict for molsimplify to structure folder
        shutil.move(self.args.name + '.dict', out_dir)

    def run_search(self):
        """
        Collate the atom and identifier data
        :return:
        """
        # Start time
        start = time.time()
        # Load pharmacophore query
        print('Loading Query')
        model = Pharmacophore.Query.from_file(self.args.catalophore)
        # Load features used in the query
        features = model.features
        catsd_features = []
        # Find all features derived from a catsd feature point and create a list containing those features
        for feature in features:
            if 'catsd' in feature.identifier:
                catsd_features.append(feature)
            else:
                pass
        # Initiate searcher
        searcher = Pharmacophore.Search()
        print('Applying Settings')
        # Apply settings
        searcher.settings.n_threads = self.args.threads
        searcher.settings.max_hit_structures = self.args.max_hit_structures
        searcher.settings.max_hits_per_structure = 1
        searcher.settings.max_hit_rmsd = self.args.rmsd
        searcher.settings.three_cubed_packing = True
        searcher.settings.complete_small_molecules = True
        # Apply filters/annotations
        print('Adding Annotations')
        model.add_feature(Pharmacophore.AnnotationFilter('is_organic', 'True'))
        # Do search
        print('Conducting Search')
        hits = searcher.search(model, self.database)
        # Obtain data for csv file and write mol file per hit
        print('Extracting Data')
        # Initiate a csd reader instance
        csd_reader = EntryReader('CSD')
        # Write csv of SMILES
        with open('SMILES' + '.csv', 'w+') as of:
            csv_writer = csv.writer(of)
            csv_writer.writerow(['Identifier', 'SMILES'])
        # Create lists for hit treatment
        list_of_data = []
        not_elements = ['Br', 'Cl', 'I', 'Li', 'Na', 'K', 'Ca', 'Mg', 'Be']
        # Stop if no hits
        if len(hits) == 0:
            print('No hits found')
            quit()
        hit_smiles = []
        failed = []
        print('Total hits: ' + str(len(hits)))
        # For each representative cluster hit perform analysis
        for h in tqdm(hits):
            # Transform hit to csd molecule instance
            mol = h.molecule
            # Molecular weight filter
            if mol.molecular_weight >= self.args.max_molecular_weight:
                continue
            else:
                pass
            # Ignore hits containing specific elements
            invalid_elements = False
            for atom in mol.atoms:
                if atom.atomic_symbol in not_elements:
                    invalid_elements = True
                    continue
                else:
                    pass
            if invalid_elements:
                continue
            # Make a copy and edit/add hydrogens
            edit_mol = mol.copy()
            edit_mol.assign_bond_types(which='unknown')
            edit_mol.add_hydrogens(mode='missing')
            edit_mol.set_formal_charges()
            # Skip if hit is a duplicate structure
            if edit_mol.smiles in hit_smiles:
                continue
            else:
                pass
            hit_smiles.append(edit_mol.smiles)
            data = [h.identifier, edit_mol.smiles]
            # Write SMILES to file
            with open('SMILES' + '.csv', 'a+') as of:
                csv_writer = csv.writer(of)
                csv_writer.writerow(data)
            # Define filename
            filename = str(h.identifier) + '.mol'
            # Get CSD entry data
            try:
                mol_entry = csd_reader.entry(mol.identifier[:mol.identifier.rindex('_')])
            except:
                mol_entry = csd_reader.entry(mol.identifier)
            # Get molecule
            csd_mol = mol_entry.molecule
            # For all catsd features find the atom index associated with the feature point
            atm_idxs = []
            for feature in catsd_features:
                # Get atom points
                atom_point = h.feature_points(feature)
                # Get xyz coordinates of the feature point
                point_coords = get_point_coords(atom_point)
                # Match atom coordinates to point feature coordinates
                atm_idx = self.get_atom_index(csd_mol, point_coords)
                atm_idxs.append(atm_idx)
            # Search hit mol for atom with same index
            coord_idxs = []
            if None in atm_idxs:
                failed.append(h.identifier)
                continue
            else:
                try:
                    for atom in atm_idxs:
                        spec_atom = edit_mol.atom(atom)
                        coord_idxs.append(spec_atom.index)
                except RuntimeError:
                    failed.append(h.identifier)
                    continue
            # Write mol2 file with structure
            with MoleculeWriter(filename, format='mol') as writer:
                writer.write(edit_mol)
            # Create data list for csv file and molsimplify ligand dictionary
            try:
                struc_data = [h.identifier[:h.identifier.rindex('_')],
                              h.identifier[h.identifier.rindex('_'):].replace('_', ''), mol_entry.chemical_name,
                              filename, str(coord_idxs), '1', h.overlay_rmsd]
                dict_data = [h.identifier[:h.identifier.rindex('_')]+':' + filename, h.identifier,
                             str(coord_idxs[0]) + ' ' + str(coord_idxs[1]), 'build custom custom', 'BA',
                             edit_mol.formal_charge]
            except:
                struc_data = [h.identifier, '1', mol_entry.chemical_name, filename, str(coord_idxs), '1',
                              h.overlay_rmsd]
                dict_data = [h.identifier+':' + filename, h.identifier+'_1', str(coord_idxs[0]) + ' ' +
                             str(coord_idxs[1]), 'build custom custom', 'BA', edit_mol.formal_charge]
            list_of_data.append(struc_data)
            # Write molsimplify dict entry
            with open(self.args.name + '.dict', 'a+') as f:
                csv_writer = csv.writer(f)
                csv_writer.writerow(dict_data)

        print('Total unique hits: ' + str(len(hit_smiles)))
        print('Failed hits: ' + str(len(failed)) + str(failed))
        print('Completed in %.2f' % (time.time() - start))
        return list_of_data

    @staticmethod
    def get_atom_index(csd_mol, point_coords):
        """
        Find atom with coordinated matching the point feature
        :param csd_mol: csd entry mol object
        :param point_coords: coordinates of the point feature
        :return: atom index
        """
        # Get atom with mol types
        atoms = csd_mol.atoms
        for atom in atoms:
            # Skip if atom has no coordinates
            if atom.coordinates is None:
                pass
            # Else truncate coordinates to 1 decimal place
            else:
                atom_coords = ', '.join('%.1f' % round(Decimal(float(a)), 1) for a in atom.coordinates)
                # Convert to list of floats
                atom_coords = [float(x) for x in atom_coords.split(',')]
                # Loops through atoms and match the feature point coords to atom coords and return the atom label of
                # the matching atom
                for coords in point_coords:
                    # Match all three coordinates with the point coordinates with a 0.1 tolerance
                    if math.isclose(atom_coords[0], coords[0], rel_tol=0.1) and \
                            math.isclose(atom_coords[1], coords[1], rel_tol=0.1) and \
                            math.isclose(atom_coords[2], coords[2], rel_tol=0.1):
                        # Get atomic symbol of matched atom
                        match_atom = atom.atomic_symbol
                        # If hydrogen get the bonded atom
                        if match_atom == 'H':
                            atom_bond = str(atom.bonds[0])
                            # Get the label of the bonded atom
                            m = re.search(r"\((\w+)\)\)", atom_bond)
                            return str(m.group(1))
                        # Else return the atom label
                        else:
                            return atom.label
                else:
                    pass

        return


if __name__ == '__main__':
    runner = Runner()
    runner.run()
