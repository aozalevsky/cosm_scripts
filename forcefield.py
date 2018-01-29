class ForceField(object):

    # Basic lattice properties
    SQ = False  # Presume default lattice type is honeycomb
    HELIXDIST = 22.0
    LIM = 21 + 6
    HC = 7

    # ------------- particles definition ---------
    count = {'T1': 1, 'T2': 2, 'T3': 3, 'T4': 4, 'T5': 5,
             'T6': 6, 'H': 7, 'T7': 7, 'T': 7, 'PT': 7, 'S': 1,
             'TT': 7, 'T1T': 1, 'T2T': 2, 'T3T': 3, 'T4T': 4,
             'T5T': 5, 'T6T': 6, 'T7T': 7, 'O': 1, 'OT': 1, 'N': 1,
             'B1': 1, 'B2': 2, 'B3': 3, 'B4': 4, 'B5': 5, 'B6': 6,
             'B7': 7, 'B1T': 1, 'B2T': 2, 'B3T': 3, 'B4T': 4,
             'B5T': 5, 'B6T': 6, 'B7T': 7, 'B': 7, 'BT': 7}
    term = ['T', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'TT', 'T1T', 'T2T',
            'T3T', 'T4T', 'T5T', 'T6T', 'T7T', 'T7',
            'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7',
            'B1T', 'B2T', 'B3T', 'B4T', 'B5T', 'B6T', 'B7T', 'BT']

    def __init__(self, lattice='h'):
        self.set_lattice_params()

    def set_sq(self, sq=False):
        self.SQ = True

        self.LIM = 36 + 7
        self.HC = 8

        self.count['PT'] = 8
        self.count['T'] = 8
        self.count['B'] = 8
        self.count['BT'] = 8
        self.count['H'] = 8
        self.count['TT'] = 8

    def set_lattice_params(self, lattice='h'):
        if lattice[0] == 's':
            self.seq_sq()
