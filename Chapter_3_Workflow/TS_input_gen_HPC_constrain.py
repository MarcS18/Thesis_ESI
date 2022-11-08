from catsd.parse.Chem import *
from catsd.parse.molsimplify import charge_spin
from catsd.parse.xtb import *
from catsd.parse.orca import *

"""
Generate HPC input files from molsimplify generated structures
xTB - input files .CHRG and .UHF
ORCA - OptTS, Freq and Energy .txt input files
.sh files are commented out as jobs are run in task arrays
CCSD are commented out as they are not required
"""


def main():

    nu_types = ['pip', 'pyr']
    file = 'terachem_input'
    charge, spin_mult = charge_spin(file)
    for file in os.listdir():
        if file.endswith('.xyz'):
            name = file.replace('.xyz', '')
            for nu in nu_types:
                if nu in name:
                    nu_name = nu
                else:
                    pass
            ts_flag = if_ts(name)
            if ts_flag:
                xtb_input(charge, spin_mult)
                xtb_input_ts(name)
                #xtb_input_hpc(name, ts_flag=True)
                xyz_name = name + '.xtbopt.xyz'
                orca_input_ts(name, xyz_name, charge, spin_mult)
                xyz_name = name + '_orca_optts.xyz'
                orca_input_freq(name, xyz_name, charge, spin_mult)
                orca_input_energy(name, xyz_name, charge, spin_mult)
                #orca_input_ccsd(name, xyz_name, charge, spin_mult)
            if not ts_flag:
                xtb_input(charge, spin_mult)
                #xtb_input_hpc(name)
                xyz_name = name + '.xtbopt.xyz'
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


def xtb_input_hpc(name, ts_flag=False):
    """
    Write sh file for running xtb calculation on ARC4
    :param name: name of the file
    :return:
    """
    file_name = name + '_xtb'
    if ts_flag:
        cmd = 'xtb --input \"constraints.inp\" \"' + name + ".xyz\"" + " --opt vtight --alpb dmf --namespace \"" + name + "\" > \"" + name + ".out\""
    else:
        cmd = 'xtb \"' + name + ".xyz\"" + " --opt vtight --alpb dmf --namespace \"" + name + "\" > \"" + name + ".out\""
    with open(file_name+'.sh', 'w+') as f:
        print('#$ -cwd -V\n'
              '#$ -l h_vmem=' + '1G' + '\n'
              '#$ -l h_rt=' + '00:30:00' + '\n'
              '#$ -l disk=' + '1G' + '\n'
              '#$ -pe smp ' + '4' + '\n'
              '#$ -m be\n' +
              cmd + '\n', file=f)


def xtb_input(charge, spin_mult):
    """
    Write property files for optimization with xtb
    :return: None
    """
    # Get charge and spin data from molsimplify generated file
    un_elec = spinmult_to_elec(spin_mult)  # Convert spin multiplicity to number of unpaired electrons
    charge_input(charge)
    spin_input(un_elec)

    return None


def xtb_input_ts(name):
    """
    Write constraints for constrained transition state optimization.
    :param name: Transition state names
    :return:
    """
    if 'TSOA' in name:
        constraints = "distance: 1, 5, auto\ndistance: 5, 14, auto\ndistance: 1, 14, auto\nangle: 1, 5, 14, auto\n"
    elif 'TSSig' in name:
        constraints = "distance: 2, 4, auto\ndistance: 4, 3, auto\nangle: 2, 4, 3, auto\n"
    else:
        return "Unsupported transition state"
    ts_constraints_constrain(constraints)

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
        print('*xyzfile '+charge+' '+spin+' '+xyz_name+'\n', file=f)
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
        print('*xyzfile '+charge+' '+spin+' '+xyz_name+'\n', file=f)
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
        print('*xyzfile ' + charge + ' ' + spin + ' ' + xyz_name + '\n', file=f)
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
        print('*xyzfile ' + charge + ' ' + spin + ' ' + xyz_name + '\n', file=f)
    #write_sh(file_name)
    return None


if __name__ == '__main__':
    main()
