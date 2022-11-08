import cclib


class Orca:
    @staticmethod
    def homo(data):
        """
        Energy of the HOMO
        :param data: ORCA output file object from cclib.ccread
        :return: Energy of the HOMO in eV
        """
        try:
            # Number of occupied orbitals
            homos = data.homos
            molecular_orbital_energies = data.moenergies
            # Energy of the HOMO
            return molecular_orbital_energies[0][homos[0]]
        except:
            return 'No Data'

    @staticmethod
    def lumo(data):
        """
        Energy of the LUMO
        :param data: ORCA output file object from cclib.ccread
        :return: Energy of the LUMO in eV
        """
        try:
            # Number of occupied orbitals
            homos = data.homos
            molecular_orbital_energies = data.moenergies
            # Energy of the HOMO
            return molecular_orbital_energies[0][homos[0]+1]
        except:
            return 'No Data'

    @staticmethod
    def imaginary_frequency(freq_data):
        """
        Return the smallest imaginary frequency
        :param freq_data: ORCA frequency output data as a cclib.ccread object
        :return: value of the imaginary frequency
        """
        try:
            return freq_data.vibfreqs[0]
        except:
            return 'No Data'

    @staticmethod
    def mulliken_charge(data, atom):
        """
        Return to muliken charge for the atom
        :param data: ORCA energy output data as a cclib.ccread object
        :param atom: Atom to get the charge for
        :return: mulliken charge of the atom
        """
        try:
            atom_charges = data.atomcharges
            # mulliken charge of the first atom in the structure
            return atom_charges['mulliken'][atom]
        except:
            return 'No Data'

    @staticmethod
    def lowdin_charge(data, atom):
        """
        Get Lowdin charge of the atom
        :param data: ORCA energy output as a cclib.ccread object
        :param atom: Atom to get the charge for
        :return: Lowdin charge of the atom
        """
        try:
            atom_charges = data.atomcharges
            # mulliken charge of the first atom in the structure
            return atom_charges['lowdin'][atom]
        except:
            return 'No Data'