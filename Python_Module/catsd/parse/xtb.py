"""
xtb.py
Functions for parsing xtb output files and retrieving stored data.
"""


def get_energy(outfile):
    """
    Return the total energy of the structure from a successful xtb calculation
    :param outfile: xtb output file
    :return: Total Energy
    """
    with open(outfile, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'total E' in line:
                return float(line.split()[-1])
            if 'TOTAL ENERGY' in line:
                return float(line.split()[-3])
    return None


def optimization_converged(outfile):
    """
    Return if geometry optimization has converged
    :param outfile: xtb output filename
    :return: True or False
    """
    with open(outfile, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'GEOMETRY OPTIMIZATION CONVERGED' in line:
                return True
    return False


def ts_constraints_fix(atoms):
    """
    Write constraint file for the constrained optimization of the ts structures generated from molsimplify
    :param atoms: Atoms to be constrained (first set including metal and substrates). Fixed for a set of substrates.
    Fixes all atom substrates. Give a higher energy structure.
    :return: None
    """
    with open("constraints.inp", "w+") as f:
        f.write('$fix\n')
        f.write('atoms: '+str(atoms)+'\n')
        f.write('$end\n')
    return None


def ts_constraints_constrain(constraints):
    """
    Write constraint file for the constrained optimization of the ts structures generated from molsimplify
    :param constraints: Constraints to be used to a transition state (first set including metal and substrates). Fixed for a TS.
    Allows other atoms to move in the substrates while fixing the TS mode. Gives lower energy structure.
    :return: None
    """
    with open("constraints.inp", "w+") as f:
        f.write('$constrain\n')
        f.write(constraints)
        f.write('$end\n')
    return None


def charge_input(charge):
    """
    Write hidden file containing charge
    :param charge: charge of complex
    :return: None
    """
    with open(".CHRG", 'w+') as of:
        of.write('%.0d\n' % (int(charge)))
    return None


def spin_input(spin):
    """
    Write hidden file containing spin
    :param spin: Number of unpaired electrons
    :return: None
    """
    with open(".UHF", 'w+') as of:
        of.write('%.0d\n' % int(spin))
    return None


def write_sh_xtb(name, memory='1G', time='00:30:00', disk='1G', nproc=4, path='xtb'):
    with open(name + '.sh', 'w+') as f:
        print('#$ -cwd -V\n'
              '#$ -l h_vmem=' + str(memory) + '\n'
              '#$ -l h_rt=' + str(time) + '\n'
              '#$ -l disk=' + str(disk) + '\n'
              '#$ -pe smp ' + str(nproc) + '\n'
              '#$ -m be\n' +
              path + ' \"' + name + '.txt\" > \"' + name + '.out\"\n', file=f)
    return None


def xtb_calc_time(outfile):
    """
    Retrieve total cpu calculation time (single core) from an output file in the format (HH:MM:SS)
    :param outfile: Calculation output file
    :return: Calculation time (single core) in HH:MM:SS else 'No Value if failed calculation
    """
    with open(outfile, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'cpu-time' in line:
                time_data = line.split(' ')
                time_dat = []
                for time in time_data:
                    try:
                        float(time)
                        time_dat.append(time)
                        time_string = str(str(int(time_dat[0] * 24) + int(time_dat[1])) + ':' + time_dat[2] + ':' + time_dat[3])
                    except:
                        pass
                return time_string
            else:
                pass

    return 'No Value'
