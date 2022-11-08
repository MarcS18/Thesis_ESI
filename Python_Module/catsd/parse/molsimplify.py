
bool_check = ('True', 'true', 'False', 'false')


def rundir(dir, inp_file):
    """Print run directory"""
    print('-rundir '+str(dir), file=inp_file)
    return


def complex_name(name, inp_file):
    """Print name of the generated complex"""
    print('-name '+str(name), file=inp_file)
    return


def core(complex_core, inp_file):
    """Print core of the complex e.g. copper"""
    print('-core '+str(complex_core), file=inp_file)
    return


def oxstate(oxidation_state, inp_file):
    """Print the oxidation state of the metal"""
    print('-oxstate '+str(oxidation_state), file=inp_file)
    return


def coord(coord_geo, inp_file):
    """Print the coordination geometry of the complex"""
    print('-coord '+str(coord_geo), file=inp_file)
    return


def geometry(geo, inp_file):
    """Print the geometry of the complex"""
    print('-geometry '+str(geo), file=inp_file)
    return


def lig(ligands, inp_file):
    """Print the ligands to be inserted into the complex"""
    print('-lig '+str(ligands), file=inp_file)
    return


def lig_freq(ligand_frequency, inp_file):
    """Print the frequency of the added ligands"""
    print('-ligocc '+str(ligand_frequency), file=inp_file)
    return


def coord_atoms(smicat, inp_file):
    """Print the coordinaing atoms"""
    print('-smicat '+str(smicat), file=inp_file)
    return


def spin(comp_spin, inp_file):
    """Print the spin of the complex"""
    print('-spin '+str(comp_spin), file=inp_file)
    return


def ff(force_field, inp_file):
    """Print the force field to be used during generation"""
    if force_field in ('uff', 'mmff94', 'gchemical', 'gaff'):
        print('-ff '+str(force_field), file=inp_file)
    else:
        return 'Invalid force field'
    return


def ffoption(force_field_option, inp_file):
    """Print when to optimize with the force field b=before, a=after, ba=before and after"""
    if force_field_option in ('b', 'a', 'ba', 'c', 'bc', 'before', 'after', 'no'):
        print('-ffoption '+str(force_field_option), file=inp_file)
    else:
        return 'Invalid force field implementation'
    return


def keepHs(keepH, inp_file):
    """Print how to handle hydrogen atoms upon coordination"""
    print('-keepHs '+str(keepH), file=inp_file)
    return


def ligalign(lig_al, inpfile):
    """Print whether to employ ligand alignment"""
    if lig_al in bool_check:
        print('-ligalign '+str(lig_al), file=inpfile)
    else:
        return 'Incorrect value for ligalign'
    return


def skipML(skipANN, inp_file):
    """Print whether to skin machine learning bond distance prediction (True or False)"""
    if skipANN in bool_check:
        print('-skipANN '+str(skipANN), file=inp_file)
    else:
        return 'Incorrect value for skipML'
    return


def replig(replace_ligand, inp_file):
    """Print whether to replace a ligand in the core"""
    if replace_ligand in bool_check:
        print('-replig '+str(replace_ligand), file=inp_file)
    else:
        return 'Incorrect value for replig'
    return


def ccatoms(coord_atoms, inp_file):
    """Print the atoms to begin replacement on the core"""
    print('-ccatoms '+str(coord_atoms), file=inp_file)
    return


def ligloc(ligand_location, inp_file):
    """Print whether to enforce ligand location criteria"""
    if ligand_location in bool_check:
        print('-ligloc '+str(ligand_location), file=inp_file)
    else:
        return 'Incorrect value for ligloc'
    return


def charge_spin(inp_file):
    """Take generated terachem input file and extract charge and spin"""
    if inp_file == 'terachem_input':
        with open(inp_file, 'r') as f:
            for line in f:
                if 'charge' in line:
                    charge_line = line.split()
                    charge = str(charge_line[1])
                elif "spinmult" in line:
                    spin_line = line.split()
                    complex_spin = str(spin_line[1])
                else:
                    pass
    else:
        return 'Incorrect file name'

    return charge, complex_spin
