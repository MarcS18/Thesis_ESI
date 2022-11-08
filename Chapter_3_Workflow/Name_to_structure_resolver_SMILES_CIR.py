import cirpy as cir
import pandas as pd
import os
import shutil
import glob


def main():

    # Import ligand from csv file from remove duplicates
    csvfile = input("Enter file name")
    smiles_cir(csvfile)


def smiles_cir(csvfile):

    lig_directory = str('Ligand_SMILES_files')
    os.makedirs(str('Ligand_SMILES_files'))
    infile = pd.read_csv(csvfile)
    ligand_column = infile['Ligand(s)'].drop_duplicates(keep="first")

    lig_failed = []

    for ligand in ligand_column:
        if ligand == 'None':
            pass
        else:
            try:
                structure = cir.resolve(ligand, 'smiles')
                with open(ligand+'.smi', 'w+') as f:
                    f.write(structure)
                    f.close()
            except:
                lig_failed.append(ligand)
                with open('Failed_Searches'+'.txt', 'w+') as f:
                    for item in lig_failed:
                        f.write("%s\n" % item)
                    f.close()

    input_files = glob.glob(str('*.smi'))
    for file in input_files:
        shutil.move(file, lig_directory)
    shutil.move('Failed_Searches.txt', lig_directory)


if __name__ == '__main__':
    main()
