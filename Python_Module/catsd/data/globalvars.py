# Dictionary containing atomic mass, atomic number, covalent radius, num valence electrons
# Data from http://www.webelements.com/ (last accessed May 13th 2015)
amassdict = {'X': (1.0, 0, 0.77, 0),     'H': (1.0079, 1, 0.37, 1),     'He': (4.002602, 2, 0.46, 2),
             'Li': (6.94, 3, 1.33, 1),   'Be': (9.0121831, 4, 1.02, 2), 'B': (10.83, 5, 0.85, 3),
             'C': (12.0107, 6, 0.77, 4), 'N': (14.0067, 7, 0.75, 5),    'O': (15.9994, 8, 0.73, 6),
             'F': (18.9984, 9, 0.71, 7), 'Ne': (20.1797, 10, 0.67, 8),  'Na': (22.99, 11, 1.55, 1),
             'Mg': (24.30, 12, 1.39, 2), 'Al': (26.98, 13, 1.26, 3),    'Si': (28.08, 14, 1.16, 4),
             'P': (30.9738, 15, 1.06, 5),'S': (32.065, 16, 1.02, 6),    'Cl': (35.453, 17, 0.99, 7),
             'Ar': (39.948, 18, 0.96, 8),'K': (39.10, 19, 1.96, 1),     'Ca': (40.08, 20, 1.71, 2),
             'Sc': (44.96, 21, 1.7, 3),  'Ti': (47.867, 22, 1.36, 4),   'V': (50.94, 23, 1.22, 5),
             'Cr': (51.9961, 24, 1.27, 6),'Mn': (54.938, 25, 1.39, 7),  'Fe': (55.84526, 26, 1.25, 8),
             'Co': (58.9332, 27, 1.26, 9), 'Ni': (58.4934, 28, 1.21, 10),'Cu': (63.546, 29, 1.38, 11),
             'Zn': (65.39, 30, 1.31, 12), 'Ga': (69.72, 31, 1.24, 3), 'Ge': (72.63, 32, 1.21, 4),
             'As': (74.92, 33, 1.21, 5),'Se': (78.96, 34, 1.16, 6),'Br': (79.904, 35, 1.14, 7),
             'Kr': (83.798, 36, 1.17, 8),'Rb': (85.47, 37, 2.10, 1), 'Sr': (87.62, 38, 1.85, 2),
             'Y': (88.91, 39, 1.63, 3), 'Zr': (91.22, 40, 1.54, 4), 'Nb': (92.91, 41, 1.47, 5),
             'Mo': (95.96, 42, 1.38, 6),'Tc': (98.9, 43, 1.56, 7), 'Ru': (101.1, 44, 1.25, 8),
             'Rh': (102.9, 45, 1.25, 9),'Pd': (106.4, 46, 1.20, 10),'Ag': (107.9, 47, 1.28, 11),
             'Cd': (112.4, 48, 1.48, 12),'In': (111.818, 49, 1.42, 3), 'Sn': (118.710, 50, 1.40, 4),
             'Sb': (121.760, 51, 1.40, 5),'Te': (127.60, 52, 1.99, 6),'I': (126.90447, 53, 1.40, 7),
             'Xe': (131.293, 54, 1.31, 8),'Cs': (132.9055, 55, 2.32, 1), 'Ba': (137.327, 56, 1.96, 2),
             'La': (138.9, 57, 1.69, 3), 'Ce': (140.116, 58, 1.63, 4), 'Pr': (140.90766, 59, 1.76, 5),
             'Nd': (144.242, 60, 1.74, 6),'Pm': (145, 61, 1.73, 7), 'Sm': (150.36, 62, 1.72, 8),
             'Eu': (151.964, 63, 1.68, 9),'Gd': (157.25, 64, 1.69, 10),'Tb': (158.92535, 65, 1.68, 11),
             'Dy': (162.500, 66, 1.67, 12), 'Ho': (164.93033, 67, 1.66, 13),'Er': (167.259, 68, 1.65, 14),
             'Tm': (168.93422, 69, 1.64, 15), 'Yb': (173.045, 70, 1.70, 16),'Lu': (174.9668, 71, 1.62, 3),
             'Hf': (178.5, 72, 1.50, 8), 'Ta': (180.9, 73, 1.38, 5), 'W': (183.8, 74, 1.46, 6),
             'Re': (186.2, 75, 1.59, 7),'Os': (190.2, 76, 1.28, 8), 'Ir': (192.2, 77, 1.37, 9),
             'Pt': (195.1, 78, 1.23, 10),'Au': (197.0, 79, 1.24, 11),'Hg': (200.6, 80, 1.49, 2),
             'Tl': (204.38, 81, 1.44, 3), 'Pb': (207.2, 82, 1.44, 4), 'Bi': (208.9804, 83, 1.51, 5),
             'Po': (208.98, 84, 1.90, 6), 'At': (209.99, 85, 2.00, 7), 'Rn': (222.6, 86, 142, 4),
             'Fr': (223.02, 87, 3.48, 8),'Ra': (226.03, 88, 2.01, 2), 'Ac': (277, 89, 1.86, 3),
             'Th': (232.0377, 90, 1.75, 4),'Pa': (231.04, 91,2.00, 5),'U': (238.02891, 92, 1.70, 6),
             'Np': (237.05, 93, 1.90, 7), 'Pu': (244.06, 94, 1.75, 8),'Am': (243.06, 95,1.80, 9),
             'Cm': (247.07, 96, 1.69, 10), 'Bk': (247.07, 97, 1.68, 11),'Cf': (251.08, 98, 1.68, 12)}

# Pa and onward should be checked

# van der Waals radii for elements
# Data from DOI: 10.1039/C3DT50599E, Dalton Trans., 2013, 42, 8617-8636
vdwrad = {'H': 1.2, 'He': 1.43, 'Li': 2.12, 'Be': 1.98, 'B': 1.91,
          'C': 1.77, 'N': 1.66, 'O': 1.50, 'F': 1.46, 'Ne': 1.58, 'Na': 2.50,
          'Mg': 2.51, 'Al': 2.25, 'Si': 2.19, 'P': 1.90, 'S': 1.89,
          'Cl': 1.82, 'Ar': 1.83, 'K': 2.73, 'Ca': 2.62, 'Sc': 2.58,
          'Ti': 2.46, 'V': 2.42, 'Cr': 2.45, 'Mn': 2.45, 'Fe': 2.44,
          'Co': 2.40, 'Ni': 2.40, 'Cu': 2.38, 'Zn': 2.39, 'Ga': 2.32,
          'Ge': 2.29, 'As': 1.88, 'Se': 1.82, 'Br': 1.86, 'Kr': 2.25,
          'Rb': 3.21, 'Sr': 2.84, 'Y': 2.75, 'Zr': 2.52, 'Nb': 2.56,
          'Mo': 2.45, 'Tc': 2.44, 'Ru': 2.46, 'Rh': 2.44, 'Pd': 2.15,
          'Ag': 2.53, 'Cd': 2.49, 'In': 2.43, 'Sn': 2.42, 'Sb': 2.47,
          'Te': 1.99, 'I': 2.04, 'Xe': 2.06, 'Cs': 3.48, 'Ba': 3.03,
          'La': 2.98, 'Ce': 2.88, 'Pr': 2.92, 'Nd': 2.95, 'Sm': 2.90,
          'Eu': 2.87, 'Gd': 2.83, 'Tb': 2.79, 'Dy': 2.87, 'Ho': 2.81,
          'Er': 2.83, 'Tm': 2.79, 'Yb': 2.80, 'Lu': 2.74, 'Hf': 2.63,
          'Ta': 2.53, 'W': 2.57, 'Re': 2.49, 'Os': 2.48, 'Ir': 2.41,
          'Pt': 2.29, 'Au': 2.32, 'Hg': 2.45, 'Tl': 2.47, 'Pb': 2.60,
          'Bi': 2.54, 'Ac': 2.8, 'Th': 2.93, 'Pa': 2.88, 'U': 2.71,
          'Np': 2.82, 'Pu': 2.81, 'Am': 2.83, 'Cm': 3.05, 'Bk': 3.4,
          'Cf': 3.05, 'Es': 2.7 }

# Period definitions for all element symbols
# Data from https://en.wikipedia.org/wiki/Group_(periodic_table) (last accessed Sept. 12th 2019)

period_1 = ['H', 'He']

period_2 = ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']

period_3 = ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']

period_4 = ['K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
            'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']

period_5 = ['Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
            'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sn', ' Te', 'I', 'Xe']

period_6 = ['Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
            'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re',
            'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']

period_7 = ['Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
            'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',
            'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

periods_dict = {'period_1': period_1, 'period_2': period_2, 'period_3': period_3,
                'period_4': period_4, 'period_5': period_5, 'period_6': period_6,
                'period_7': period_7}

# Group definitions for all element symbols
# Data from https://en.wikipedia.org/wiki/Group_(periodic_table) (last accessed Sept. 12th 2019)

hydrogen = ['H']  # Note H not typically included in group 1

group_1 = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']

group_2 = ['Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra']

group_3 = ['Sc', 'Y']  # Some IUPAC Initiatives to call either 'La' and 'Ac' grp 3 or 'Lu' and 'Lr

group_4 = ['Ti', 'Zr', 'Hf', 'Rf']

group_5 = ['V', 'Nb', 'Ta', 'Db']

group_6 = ['Cr', 'Mo', 'W', 'Sg']

group_7 = ['Mn', 'Tc', 'Re', 'Bh']

group_8 = ['Fe', 'Ru', 'Os', 'Hs']

group_9 = ['Co', 'Rh', 'Ir', 'Mt']

group_10 = ['Ni', 'Pd', 'Pt', 'Ds']

group_11 = ['Cu', 'Ag', 'Au', 'Rg']

group_12 = ['Zn', 'Cd', 'Hg', 'Cn']

group_13 = ['B', 'Al', 'Ga', 'In', 'Tl', 'Nh']

group_14 = ['C', 'Si', 'Ge', 'Sn', 'Pb', 'Fl']

group_15 = ['N', 'P', 'As', 'Sb', 'Bi', 'Mc']

group_16 = ['O', 'S', 'Se', 'Te', 'Po', 'Lv']

group_17 = ['F', 'Cl', 'Br', 'I', 'At', 'Ts']

group_18 = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn', 'Og']

lanthanides = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
               'Ho', 'Er', 'Tm', 'Yb', 'Lu']

actinides = ['Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
             'Es', 'Fm', 'Md', 'No', 'Lr']

halides = {
    'iodide': {'name': 'iodide',
               'smiles': '[I-]',
               'catom': '[1]'},
    'bromide': {'name': 'bromide',
                'smiles': '[Br-]',
                'catom': '[1]'},
    'chloride': {'name': 'chloride',
                 'smiles': '[Cl-]',
                 'catom': '[1]'},
}

# Dictionary of nucleophiles
nucleophiles = {
    'pip': {'name': 'piperidine',
            'smiles': 'C1CCNCC1',
            'catom': '[4]'},
    'pyr': {'name': 'pyrrolidinone',
            'smiles': 'C1CCC(N1)=O',
            'catom': '[5]'}
}

ts_cores = {
    'TSOA_pip': {'name': 'TSOA_pip',
                 'core_name': 'tsoa_pip',
                 'ccatoms': '6,15'},
    'TSSig_pip': {'name': 'TSOA_pip',
                  'core_name': 'tssig_pip',
                  'ccatoms': '4,15'},
    'TSOA_pyr': {'name': 'TSOA_pyr',
                 'core_name': 'tsoa_pyr',
                 'ccatoms': '6,15'},
    'TSSig_pyr': {'name': 'TSSig_pyr',
                  'core_name': 'tssig_pyr',
                  'ccatoms': '4,15'}
}


groups_dict = {'group_1': group_1, 'group_2': group_2, 'group_3': group_3,
               'group_4': group_4, 'group_5': group_5, 'group_6': group_6,
               'group_7': group_7, 'group_8': group_8, 'group_9': group_9,
               'group_10': group_10, 'group_11': group_11, 'group_12': group_12,
               'group_13': group_13, 'group_14': group_14, 'group_15': group_15,
               'group_16': group_16, 'group_17': group_17, 'group_18': group_18,
               'lanthanides': lanthanides, 'actinides': actinides, 'hydrogen': hydrogen}

# Metals (includes alkali, alkaline earth, and transition metals)
alkali_and_alkaline_earth = ['Li', 'li', 'LI', 'lithium', 'Be', 'be', 'BE', 'beryllium',
                             'Na', 'na', 'NA', 'sodium', 'Mg', 'mg', 'MG', 'magnesium',
                             'Al', 'al', 'AL', 'aluminum', 'aluminium',
                             'K', 'k', 'potassium', 'Ca', 'ca', 'CA', 'calcium',
                             'Rb', 'rb', 'RB', 'rubidium', 'Sr', 'sr', 'SR', 'strontium',
                             'Cs', 'cs', 'CS', 'cesium', 'Ba', 'ba', 'BA', 'barium',
                             'Fr', 'fr', 'FR', 'francium', 'Ra', 'ra', 'RA', 'radium']

heavy_metals_and_metalloids = ['Ga', 'ga', 'GA', 'gallium',
                                'In', 'in', 'IN', 'indium', 'Sn', 'sn', 'SN', 'tin',
                                'Tl', 'tl', 'TL', 'thallium', 'Pb', 'pb', 'PB', 'lead',
                                'Bi', 'bi', 'BI', 'bismuth', 'Po', 'po', 'PO', 'polonium',
                                'La', 'la', 'LA', 'lanthanum',
                                'Ce', 'ce', 'CE', 'cerium', 'Pr', 'pr', 'PR', 'praseodymium',
                                'Nd', 'nd', 'ND', 'neodymium', 'Pm', 'pm', 'PM', 'promethium',
                                'Sm', 'sm', 'SM', 'samarium', 'Eu', 'eu', 'EU', 'europium',
                                'Gd', 'gd', 'GD', 'gadolinium', 'Tb', 'tb', 'TB', 'terbium',
                                'Dy', 'dy', 'DY', 'dysprosium', 'Ho', 'ho', 'HO', 'holmium',
                                'Er', 'er', 'ER', 'erbium', 'Tm', 'tm', 'TM', 'thulium',
                                'Yb', 'yb', 'YB', 'ytterbium', 'Lu', 'lu', 'LU', 'lutetium',
                                'Ac', 'ac', 'AC', 'actinium', 'Th', 'th', 'TH', 'thorium',
                                'Pa', 'pa', 'PA', 'proactinium', 'U', 'u', 'uranium',
                                'Np', 'np', 'NP', 'neptunium', 'Pu', 'pu', 'PU', 'plutonium',
                                'Am', 'am', 'AM', 'americium', 'Cu', 'cu', 'CU', 'curium',
                                'Bk', 'bk', 'BK', 'berkelium', 'Cf', 'cf', 'CF', 'californium',
                                'Es', 'es', 'ES', 'einsteinium', 'Fm', 'fm', 'FM', 'fermium',
                                'Md', 'md', 'MD', 'mendelevium', 'No', 'no', 'NO', 'nobelium',
                                'Lr', 'lr', 'LR', 'lawrencium']
metalslist = ['Li', 'li', 'LI', 'lithium', 'Be', 'be', 'BE', 'beryllium',
              'Na', 'na', 'NA', 'sodium', 'Mg', 'mg', 'MG', 'magnesium',
              'Al', 'al', 'AL', 'aluminum', 'aluminium',
              'K', 'k', 'potassium', 'Ca', 'ca', 'CA', 'calcium',
              'Rb', 'rb', 'RB', 'rubidium', 'Sr', 'sr', 'SR', 'strontium',
              'Cs', 'cs', 'CS', 'cesium', 'Ba', 'ba', 'BA', 'barium',
              'Fr', 'fr', 'FR', 'francium', 'Ra', 'ra', 'RA', 'radium',
              'Sc', 'sc', 'SC', 'scandium', 'Ti', 'ti', 'TI', 'titanium',
              'V', 'v', 'vanadium', 'Cr', 'cr', 'CR', 'chromium',
              'Mn', 'mn', 'MN', 'manganese', 'Fe', 'fe', 'FE', 'iron',
              'Co', 'co', 'CO', 'cobalt', 'Ni', 'ni', 'NI', 'nickel',
              'Cu', 'cu', 'CU', 'copper', 'Zn', 'zn', 'ZN', 'zinc',
              'Ga', 'ga', 'GA', 'gallium',
              'Y', 'y', 'yttrium', 'Zr', 'zr', 'ZR', 'zirconium',
              'Nb', 'nb', 'NB', 'niobium', 'Mo', 'mo', 'MO', 'molybdenum',
              'Tc', 'tc', 'TC', 'technetium', 'Ru', 'ru', 'RU', 'ruthenium',
              'Rh', 'rh', 'RH', 'rhodium', 'Pd', 'pd', 'PD', 'palladium',
              'Ag', 'ag', 'AG', 'silver', 'Cd', 'cd', 'CD', 'cadmium',
              'In', 'in', 'IN', 'indium', 'Sn', 'sn', 'SN', 'tin',
              'Hf', 'hf', 'HF', 'hafnium', 'Ta', 'ta', 'TA', 'tantalum',
              'W', 'w', 'tungsten', 'Re', 're', 'RE', 'rhenium',
              'Os', 'os', 'OS', 'osmium', 'Ir', 'ir', 'IR', 'iridium',
              'Pt', 'pt', 'PT', 'platinum', 'Au', 'au', 'AU', 'gold',
              'Hg', 'hg', 'HG', 'mercury', 'X',
              'Tl', 'tl', 'TL', 'thallium', 'Pb', 'pb', 'PB', 'lead',
              'Bi', 'bi', 'BI', 'bismuth', 'Po', 'po', 'PO', 'polonium',
              'La', 'la', 'LA', 'lanthanum',
              'Ce', 'ce', 'CE', 'cerium', 'Pr', 'pr', 'PR', 'praseodymium',
              'Nd', 'nd', 'ND', 'neodymium', 'Pm', 'pm', 'PM', 'promethium',
              'Sm', 'sm', 'SM', 'samarium', 'Eu', 'eu', 'EU', 'europium',
              'Gd', 'gd', 'GD', 'gadolinium', 'Tb', 'tb', 'TB', 'terbium',
              'Dy', 'dy', 'DY', 'dysprosium', 'Ho', 'ho', 'HO', 'holmium',
              'Er', 'er', 'ER', 'erbium', 'Tm', 'tm', 'TM', 'thulium',
              'Yb', 'yb', 'YB', 'ytterbium', 'Lu', 'lu', 'LU', 'lutetium',
              'Ac', 'ac', 'AC', 'actinium', 'Th', 'th', 'TH', 'thorium',
              'Pa', 'pa', 'PA', 'proactinium', 'U', 'u', 'uranium',
              'Np', 'np', 'NP', 'neptunium', 'Pu', 'pu', 'PU', 'plutonium',
              'Am', 'am', 'AM', 'americium', 'Cu', 'cu', 'CU', 'curium',
              'Bk', 'bk', 'BK', 'berkelium', 'Cf', 'cf', 'CF', 'californium',
              'Es', 'es', 'ES', 'einsteinium', 'Fm', 'fm', 'FM', 'fermium',
              'Md', 'md', 'MD', 'mendelevium', 'No', 'no', 'NO', 'nobelium',
              'Lr', 'lr', 'LR', 'lawrencium']

metals_conv = {'scandium': 'Sc', 'titanium': 'Ti', 'vanadium': 'V', 'chromium': 'Cr', 'manganese': 'Mn',
               'iron': 'Fe', 'cobalt': 'Co', 'nickel': 'Ni', 'copper': 'Cu', 'zinc': 'Zn',
               'yttrium': 'Y', 'zirconium': 'Zr', 'niobium': 'Nb', 'molybdenum': 'Mo', 'technetium': 'Tc',
               'ruthenium': 'Ru', 'rhodium': 'Rh', 'palladium': 'Pd', 'silver': 'Ag', 'cadmium': 'Cd',
               'lanthanum': 'La', 'hafnium': 'Hf',
               'tantalum': 'Ta', 'tungsten': 'W',
               'rhenium': 'Re', 'osmium': 'Os',
               'iridium': 'Ir', 'platinum': 'Pt',
               'gold': 'Au', 'mercury': 'Hg'}

# d-electron counts of transition metals
mtlsdlist = {'sc': 1, 'ti': 2, 'v': 3, 'cr': 4, 'mn': 5, 'fe': 6, 'co': 7, 'ni': 8, 'cu': 9, 'zn': 10,
             'y': 1, 'zr': 2, 'nb': 3, 'mo': 4, 'tc': 5, 'ru': 6, 'rh': 7, 'pd': 8, 'ag': 9, 'cd': 10,
             'hf': 2, 'ta': 3, 'w': 4, 're': 5, 'os': 6, 'ir': 7, 'pt': 8, 'au': 9, 'hg': 10}

# Default spins for each d-electron count (make this metal/oxidation state specific)
defaultspins = {0: '1', 1: '2', 2: '3', 3: '4', 4: '5',
                5: '6', 6: '5', 7: '4', 8: '3', 9: '2', 10: '1'}

# Elements sorted by atomic number
elementsbynum = ['H', 'He',
                 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
                 'K', 'Ca',
                 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
                 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',
                 'Xe',
                 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
                 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
                 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
                 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo']

# Electronegativity (Pauling) by atom symbol
endict = {"H": 2.20, "He": 4.16,
          "Li": 0.98, "Be": 1.57, "B": 2.04, "C": 2.55, "N": 3.04, "O": 3.44, "F": 3.98,
          "Na": 0.93, "Mg": 1.31, "Al": 1.61, "Si": 1.90, "P": 2.19, "S": 2.58, "Cl": 3.16,
          "K": 0.82, "Ca": 1.00, "Sc": 1.36, "Ti": 1.54, "V": 1.63, "Cr": 1.66,
          "Mn": 1.55, "Fe": 1.83, "Co": 1.88, "Ni": 1.91, "Cu": 1.90, "Zn": 1.65, "Ga": 1.81,
          "Ge": 2.01, "As": 2.18, "Se": 2.55, "Br": 2.96, "Rb": 0.82, "Sr": 0.95, "Y": 1.22,
          "Zr": 1.33, "Nb": 1.60, "Mo": 2.16, "Tc": 2.10, "Ru": 2.20, "Rh": 2.28,
          "Pd": 2.20, "Ag": 1.93, "Cd": 1.69, "In": 1.78, "Sn": 1.96, "Sb": 2.05, "I": 2.66,
          "Cs": 0.79, "Ba": 0.89, "Hf": 1.30, "Ta": 1.50, "W": 2.36, "Re": 1.90, "Os": 2.20, "Ir": 2.20,
          "Pt": 2.28, "Au": 2.54, "Hg": 2.00, "Tl": 1.62, "Pb": 2.33, "Bi": 2.02,
          "La": 1.10, "Ce": 1.12, "Pr": 1.13, "Nd": 1.14, "Sm": 1.17,
          "Gd": 1.20, "Dy": 1.22, "Ho": 1.23, "Er": 1.24, "Tm": 1.25, "Lu": 1.27,
          "Fr": 0.7, "Ra": 0.9, "Ac": 1.1, "Th": 1.3, "Pa": 1.5, "U": 1.38, "Np": 1.36, "Pu": 1.28,
          "Am": 1.3, "Cm": 1.3, "Bk": 1.3, "Cf": 1.3, "Es": 1.3, "Fm": 1.3, "Md": 1.3, "No": 1.3,
          "Yb": 1.1, "Eu": 1.2, "Tb": 1.1, "Te": 2.10}

# Polarizability (alpha) by atom symbol
# From https://www.tandfonline.com/doi/full/10.1080/00268976.2018.1535143
# Last accessed 4/28/20

poldict = {"H": 4.50711, "He": 1.38375,
           "Li": 164.1125, "Be": 37.74, "B": 20.5, "C": 11.3, "N": 7.4,
           "O":5.3, "F": 3.74, "Ne": 2.66, "Na": 162.7, "Mg":71.2, "Al": 57.8, "Si": 37.3, "P": 25,
           "S": 19.4, "Cl": 14.6, "Ar": 11.083, "K": 289.7, "Ca": 160.8, "Sc": 97, "Ti": 100,
           "V": 87, "Cr": 83, "Mn": 68, "Fe": 62, "Co": 55, "Ni": 49, "Cu": 46.5, "Zn": 38.67,
           "Ga": 50, "Ge": 40, "As": 30, "Se": 28.9, "Br": 21, "Kr": 16.78, "Rb": 319.8, "Sr": 197.2,
           "Y": 162, "Zr": 112, "Nb": 98, "Mo": 87, "Tc": 79, "Ru": 72, "Rh": 66, "Pd": 26.14,
           "Ag": 55, "Cd": 46, "In": 65, "Sn": 53, "Sb": 43, "Te": 38, "I": 32.9, "Xe": 27.32,
           "Cs": 400.9, "Ba": 272, "La": 215, "Ce": 205, "Pr": 216, "Nd": 208, "Pm": 200, "Sm": 192,
           "Eu": 184, "Gd": 158, "Tb": 170, "Dy": 163, "Ho": 156, "Er": 150, "Tm": 144,
           "Yb": 139, "Lu": 137, "Hf": 103, "Ta": 74, "W": 68, "Re": 62, "Os": 57, "Ir": 54,
           "Pt": 48, "Au": 36, "Hg": 33.91, "Tl": 50, "Pb": 47, "Bi": 48, "Po": 44, "At": 42,
           "Rn": 35, "Fr": 317.8, "Ra": 246, "Ac": 203, "Pa": 154, "U": 129, "Np": 151, "Pu": 132,
           "Am": 131, "Cm": 144, "Bk": 125, "Cf": 122, "Es": 118, "Fm": 113, "Md": 109, "No": 110,
           "Lr": 320, "Rf": 112, "Db": 42, "Sg": 40, "Bh": 38, "Hs": 36, "Mt": 34, "Ds": 32,
           "Rg": 32, "Cn": 28, "Nh": 29, "Fl": 31, "Mc": 71, "Ts": 76, "Og": 58}


# Roman numerals
romans = {'I': '1', 'II': '2', 'III': '3', 'IV': '4',
          'V': '5', 'VI': '6', 'VII': '7', 'VIII': '8'}

# bondsdict
bondsdict = {"H": 1, "Li": 1, "Be": 2, "B": 3, "C": 4, "N": 3, "O": 2, "F": 1,
             "Na": 1, "Mg": 2, "Al": 3, "Si": 4, "P": 3, "S": 2, "Cl": 1,
             "As": 3, "Se": 2, "Br": 1, "I": 1}

# triple bonds dictionry: Defined as 0.5*(double bond dist + triple bond dist)
# bond lengths are from http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html
tribonddict = {("C", "C"): 1.27, ("C", "N"): 1.235, ("C", "O"): 1.165, ("N", "N"): 1.175,
               ("N", "C"): 1.235, ("O", "C"): 1.165}


"""
Gibbs free energy of common reagents in DMF, xtb and B97-3c, coupled cluster and various dft methods with xtb correction
ccsd_old is the energy used in the B97-3c Energy benchmark using ORCA 4.2.1
ccsd is the new energy calculated using ORCA 5.0.2
"""

dmf_common = {
    'iodide': {'xtb': -4.145816493,
               'B97-3c': -298.2617209,
               'B97-3c_70': -298.2617209,
               'ccsd_old': -297.3170952,
               'ccsd': -297.315141,
               'r2SCAN-3c': -297.9215522,
               'PBE0': -297.9034076,
               'TPSS': -297.7776279,
               'TPSSh': -297.7994433},
    'piperidine': {'xtb': -19.12108408,
                   'B97-3c': -251.6613233,
                   'B97-3c_70': -251.6667664,
                   'ccsd_old': -251.3346061,
                   'ccsd': -251.3366574,
                   'r2SCAN-3c': -251.7052983,
                   'PBE0': -251.5831667,
                   'TPSS': -251.9339726,
                   'TPSSh': -251.9112888},
    'pyrrolidinone': {'xtb': -19.07560702,
                      'B97-3c': -286.438338,
                      'B97-3c_70': -286.4437647,
                      'ccsd_old': 0,
                      'ccsd': -286.0996381,
                      'r2SCAN-3c': -286.5060177,
                      'PBE0': -286.3609207,
                      'TPSS': -286.7364363,
                      'TPSSh': -286.7055393},
    'iodobenzene': {'xtb': -19.14473138,
                    'B97-3c': -529.5633063,
                    'B97-3c_70': -529.5691971,
                    'ccsd_old': -528.294146,
                    'ccsd': -528.2946219,
                    'r2SCAN-3c': -529.2716701,
                    'PBE0': -529.1300271,
                    'TPSS': -529.33735,
                    'TPSSh': -529.3310189},
    'N_phenyl_piperidine': {'xtb': -33.94899071,
                            'B97-3c': -482.5301061,
                            'B97-3c_70': -482.5372188,
                            'ccsd_old': -481.8758961,
                            'ccsd': -481.8807981,
                            'r2SCAN-3c': -482.6337866,
                            'PBE0': -482.384318,
                            'TPSS': -483.061158,
                            'TPSSh': -483.0122785},
    'N_phenyl_pyrrolidinone': {'xtb': -33.90011439,
                               'B97-3c': -517.3046708,
                               'B97-3c_70': -517.311698,
                               'ccsd_old': 0,
                               'ccsd': -516.6435768,
                               'r2SCAN-3c': -517.4325275,
                               'PBE0': -517.1594278,
                               'TPSS': -517.8614894,
                               'TPSSh': -517.8043933},
    'carbonate': {'xtb': -14.96153031,
                  'B97-3c': -263.9106901,
                  'B97-3c_70': -263.9152400,
                  'ccsd_old': -263.6492102,
                  'ccsd': -263.6549378,
                  'r2SCAN-3c': -263.9953859,
                  'PBE0': -263.8483653,
                  'TPSS': -264.1554823,
                  'TPSSh': -264.1216463},
    'hydrogen_carbonate': {'xtb': -15.20638829,
                           'B97-3c': -264.4434517,
                           'B97-3c_70': -264.4481104,
                           'ccsd_old': -264.1929425,
                           'ccsd': -264.196782,
                           'r2SCAN-3c': -264.5273195,
                           'PBE0': -264.3859246,
                           'TPSS': -264.6934531,
                           'TPSSh': -264.6603292},
    'caesium': {'xtb': 0.142543716,
                'B97-3c': -20.1166441,
                'B97-3c_70': -20.1166441,
                'ccsd_old': -19.94847101,
                'ccsd': -19.94339042,
                'r2SCAN-3c': -20.079266,
                'PBE0': -20.06662143,
                'TPSS': -20.07627415,
                'TPSSh': -20.07415349},
}
# Pure TPSSh/def2-TZVP Gibbs Free Energy calculated with ORCA 5.0.2 in DMF
dmf_common_tpssh = {
    'iodide': {'TPSSh': -297.79442300},
    'piperidine': {'TPSSh': -251.90998556},
    'pyrrolidinone': {'TPSSh': -286.70365767},
    'iodobenzene': {'TPSSh': -529.32921831},
    'N_phenyl_piperidine': {'TPSSh': -483.00817162},
    'N_phenyl_pyrrolidinone': {'TPSSh': -517.79885057},
    'carbonate': {'TPSSh': -264.12019921},
    'hydrogen_carbonate': {'TPSSh': -264.65981217},
    'caesium': {'TPSSh': -20.06917469},
}