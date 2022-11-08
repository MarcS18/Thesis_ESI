import datetime
import re

"""
ORCA.py
Functions for parsing ORCA output files and retrieving stored data.
"""


def get_energy(outfile):
    with open(outfile, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'FINAL SINGLE POINT ENERGY' in line:
                return float(line.split()[4])
    return 'NaN'


def get_enthalpy(outfile):
    with open(outfile, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'Total Enthalpy' in line:

                try:
                    return float(line.split()[-2])

                except ValueError:
                    break
    return 'NaN'


def get_free_energy(outfile):
    with open(outfile, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'Final Gibbs free energy' in line:
                try:
                    return float(line.split()[-2])
                except ValueError:
                    break
    return 'NaN'


def optimization_converged(outfile):
    with open(outfile, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'THE OPTIMIZATION HAS CONVERGED' in line:
                return True
    return False


def write_solvent(solvent, f):
    """
    Write solvent to ORCA input file
    :param solvent: Solvent to be used in the calculation
    :return: None
    """
    print(f'%cpcm\n'
          f'smd true\n'
          f'SMDsolvent \"{solvent}\"\n'
          f'end\n', file=f)
    return None


def write_geom(f, hess_freq=5):
    """
    Write geometry parameters to ORCA input file
    :param f: Input file name
    :return: None
    """
    print(f'%geom\n'
          f'Calc_Hess true\n'
          f'NumHess true\n'
          f'Recalc_Hess {hess_freq}\n'
          f'end\n', file=f)
    return None


def write_scf(f, maxiter=2000):
    """
    Write SCF parameters
    :param f: Input file name
    :return: None
    """
    print(f'%scf\n'
          f'MaxIter {maxiter}\n'
          f'end\n', file=f)
    return None


def write_maxcore(f, memory=8000):
    """
    Write maxcore for CCSD(T) energy calculation
    :param f: Input file name
    :param memory: Memory of the calculation
    :return:
    """
    print(f'%maxcore {int(memory * 0.75)}\n', file=f)
    return None


def write_sh(name, memory='8G', time='48:00:00', disk='1G', nproc=4, path='/apps/applications/orca/4.2.1/1/default/bin/orca'):
    with open(name + '.sh', 'w+') as f:
        print(f'#$ -cwd -V\n'
              f'#$ -l h_vmem={memory}\n'
              f'#$ -l h_rt={time}\n'
              f'#$ -l disk={disk}\n'
              f'#$ -pe smp {nproc}\n'
              f'#$ -m be\n'
              f'module add orca/4.2.1\n' +
              f"{path} \"{name}.txt\" > \"{name}.out\"\n'", file=f)
    return None


def get_correction(freq_file):
    with open(freq_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'G-E(el)' in line:
                try:
                    return float(line.split()[-4])
                except ValueError:
                    break
    return None


def orca_calc_time(outfile):
    """
    Calculate computational time as single core from orca output file
    :param outfile: Calculation output file
    :return: Time in HH:MM:SS else 'No Value' if failed calculation
    """
    with open(outfile, 'r') as f:
        lines = f.readlines()
        # Get the number of threads used in the calculation
        for line in lines:
            # Find PAL in keywords in input file
            if 'PAL' in line:
                # Get only commands line
                commands = line[:-1]
                # Split line into words
                words = commands.split(' ')
                # Get only threads keyword and extract number of threads
                for word in words:
                    if 'PAL' in word:
                        threads = word.split('L')
                        threads = threads[-1]
            else:
                pass
            # Get run time
            if 'TOTAL RUN TIME:' in line:
                time_data = line.split(' ')
                time_dat = []
                for time in time_data:
                    # Only retrieve integers in line
                    try:
                        int(time)
                        time_dat.append(time)
                        # Convert days to hours in format (HH:MM:SS)
                        time_string = str(str(int(time_dat[0] * 24) + int(time_dat[1])) + ':' + time_dat[2] + ':' + time_dat[3])
                        # Create time delta
                        my_sum = datetime.timedelta()
                        # Split wall time
                        (h, m, s) = time_string.split(':')
                        # Convert wall time to timedelta
                        d = datetime.timedelta(hours=int(h), minutes=int(m), seconds=float(s))
                        # Multiple wall time by number of threads for single core time
                        my_sum += (d * float(threads))
                        days, seconds = my_sum.days, my_sum.seconds
                        hours = days * 24 + seconds // 3600
                        minutes = (seconds % 3600) // 60
                        seconds = seconds % 60
                        my_sum = str(str(hours) + ':' + str(minutes) + ':' + str(seconds))
                    except:
                        pass
                return str(my_sum)
            else:
                pass

    return 'No Value'


def check_negative_frequencies(file, tolerance=0, ts_tolerance=-40, mode='lessthan'):
    """
    Checks for negative frequencies in an ORCA 4.2 output file against a tolerance.
    :param file: Input ORCA 4.2 .out file
    :param tolerance: Tolerance for the magnitude of the negative frequencies. Default tolerance = 0
    :param ts_tolerance: Tolerance for the magnitude of the transition state imaginary frequency. Default = -40
    :param mode: mode for checking frequencies, lessthan: smaller than tolerance, ts: greater than the tolerance (for TS's)
    :return: Minimum if no negative frequencies, Not minimum if negative frequencies
    """
    # Create empty list for file.
    data = []
    # Create empty list for frequencies
    vib_freqs = []
    with open(file, 'r') as f:
        try:
            for line in f:
                if line.strip() == 'VIBRATIONAL FREQUENCIES':
                    break
            for line in f:
                if line.strip() == 'NORMAL MODES':
                    break
                data.append(line)
            data = data[4:-3]
            for line in data:
                # Remove empty entries from split line on spaces
                split_line = list(filter(None, line.split(' ')))
                # Append wavenumber to vib_freqs
                vib_freqs.append(split_line[1])
            # Check bond for correct TS components and return correct if found
            if len(data) == 0:
                return 'Calc Fail'
            if mode == 'lessthan':
                if any(float(n) < tolerance for n in vib_freqs):
                    return 'Not Minimum'
                else:
                    return 'Minimum'
            elif mode == 'ts':
                if float(vib_freqs[6]) < ts_tolerance and float(vib_freqs[7]) > tolerance:
                    return 'Transition State'
                else:
                    return 'Not Transition State'
            else:
                return 'Unsupported Mode'
        except:
            # If not present return calculation failed
            return 'Calc Fail'
