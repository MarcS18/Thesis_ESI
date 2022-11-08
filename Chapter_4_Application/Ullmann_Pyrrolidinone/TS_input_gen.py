from catsd.parse.Chem import *
from catsd.parse.molsimplify import charge_spin
from catsd.parse.orca import *

"""
Generate HPC input files from molsimplify generated structures
ORCA - OptTS, Freq and Energy .txt input files
.sh files are commented out as jobs are run in task arrays
CCSD are commented out as they are not required
"""


def main():

    cs_file = 'terachem_input'
    charge, spin_mult = charge_spin(cs_file)
    for file in os.listdir():
        if file.endswith('.xyz'):
            name = file.replace('.xyz', '')
            ts_flag = if_ts(name)
            xyz_name = file
            if ts_flag:
                orca_input_opt(name, xyz_name, charge, spin_mult)
                xyz_name = name + '_orca_opt.xyz'
                orca_input_ts(name, xyz_name, charge, spin_mult)
                xyz_name = name + '_orca_optts.xyz'
                orca_input_freq(name, xyz_name, charge, spin_mult)
                orca_input_energy(name, xyz_name, charge, spin_mult)
                #orca_input_ccsd(name, xyz_name, charge, spin_mult)
            if not ts_flag:
                orca_input_opt(name, xyz_name, charge, spin_mult)
                xyz_name = name + '_orca_opt.xyz'
                orca_input_freq(name, xyz_name, charge, spin_mult)
                orca_input_energy(name, xyz_name, charge, spin_mult)
                #orca_input_ccsd(name, xyz_name, charge, spin_mult)
    return


def if_ts(name):
    """
    Return if structure is a transition state
    :param name: name of the file
    :return: True if transition state, False if not a transition state
    """
    ts_types = ['TSOA', 'TSSig']
    for ts in ts_types:
        if ts in name:
            return True

    return False


def orca_input_opt(name, xyz_name, charge, spin, solvent='dmf'):
    """
    Write orca preopt/opt input file
    :param name: Name of the input file
    :param xyz_name: Name of the xyz file
    :param charge: Charge of the complex
    :param spin: Spin multiplicity of the complex
    :param solvent: Solvent to be used in the calculation
    :return: None
    """
    file_name = name + '_orca_opt'
    with open(file_name+'.txt', 'w+') as f:
        print('! XTB2 PAL4 TightOpt\n', file=f)
        if 'TSOA' in name:
            print('%geom\n'
                  'MaxIter 9999\n'
                  'Constraints\n'
                  '\t{B 0 4 C}\n'
                  '\t{B 4 13 C}\n'
                  '\t{B 0 13 C}\n'
                  '\t{A 0 4 13 C}\n'
                  '\tend\n'
                  'end\n', file=f)
        elif 'TSSig' in name:
            print('%geom\n'
                  'MaxIter 9999\n'
                  'Constraints\n'
                  '\t{B 1 3 C}\n'
                  '\t{B 2 3 C}\n'
                  '\t{A 1 2 3 C}\n'
                  '\tend\n'
                  'end\n', file=f)
        else:
            print('%geom\n'
                  'MaxIter 9999\n'
                  'end\n', file=f)
        write_solvent(solvent, f)
        print(f'*xyzfile {charge} {spin} {xyz_name}\n', file=f)

    return None


def orca_input_freq(name, xyz_name, charge, spin, solvent='dmf'):
    """
    Write ORCA frequency input file
    :param name: Name of the input file
    :param xyz_name: Name of the xyz file
    :param charge: Charge of the complex
    :param spin: Spin multiplicity of the complex
    :param solvent: Solvent to be used in the calculation
    :return: None
    """
    file_name = name + '_orca_freq'
    with open(file_name+'.txt', 'w+') as f:
        print('! XTB2 PAL4 NumFreq\n', file=f)
        write_solvent(solvent, f)
        print(f'*xyzfile {charge} {spin} {xyz_name}\n', file=f)
    #write_sh(file_name, memory='1G', time='06:00:00', path='/home/home02/pmmass/orca/orca')
    return None


def orca_input_ts(name, xyz_name, charge, spin, solvent='dmf'):
    """
    Write ORCA transition state optimization file
    :param name: input file name
    :param xyz_name: name of the xyz file
    :param charge: charge of complex
    :param spin: spin multiplicity of complex
    :param solvent: solvent to be used
    :return: None
    """
    file_name = name + '_orca_optts'
    with open(file_name+'.txt', 'w+') as f:
        print('! XTB2 PAL4 OptTS\n', file=f)
        write_geom(f)
        write_solvent(solvent, f)
        print(f'*xyzfile {charge} {spin} {xyz_name}\n', file=f)
    #write_sh(file_name, memory='1G', time='24:00:00', path='/home/home02/pmmass/orca/orca')
    return None


def orca_input_energy(name, xyz_name, charge, spin, solvent='dmf'):
    """
    Generate input ORCA file for B97-3c energy calculation
    :param name: name of calculation
    :param xyz_name: name of xyz file
    :param solvent: solvent
    :param charge: charge of complex
    :param spin: spin of complex
    :return:
    """
    filename = name + '_orca_energy'
    with open(filename + '.txt', 'w+') as f:
        print('! B97-3c PAL4 TIGHTSCF SlowConv\n', file=f)
        write_scf(f)
        write_solvent(solvent, f)
        print(f'*xyzfile {charge} {spin} {xyz_name}\n', file=f)
    #write_sh(filename, memory='1G', time='24:00:00')
    return None


def orca_input_ccsd(name, xyz_name, charge, spin, solvent='dmf'):
    """
    Write DLPNO-CCSD(T) energy calculation input files, .txt and .sh for supercomputer
    :param name: name of the complex
    :param xyz_name: name of the xyz file to operate on
    :param charge: charge of the complex
    :param spin: spin of the complex
    :param solvent: solvent used in calculation
    :return:
    """
    file_name = name + '_orca_ccsd'
    with open(file_name + '.txt', 'w+') as f:
        print('! DLPNO-CCSD(T) DEF2-TZVPP DEF2-TZVPP/C PAL4 TIGHTSCF\n', file=f)
        write_maxcore(f, memory=8000)
        write_scf(f)
        write_solvent(solvent, f)
        print(f'*xyzfile {charge} {spin} {xyz_name}\n', file=f)
    #write_sh(file_name)
    return None


if __name__ == '__main__':
    main()
