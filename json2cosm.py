#!/usr/bin/env python
# encoding: utf-8

import sys
import argparse
import string
import copy
import math
import numpy as np
import traceback as tb
from collections import OrderedDict

import util
from util import COSMError
from forcefield import ForceField


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-i', '--input', required=True,
        action='store', dest='input',
        help='input cadnano json file')

    parser.add_argument(
        '-o', '--output', required=True,
        action='store', dest='output',
        help='output pdb file')

    parser.add_argument(
        '-r', '--restr', required=True,
        action='store', dest='cnct',
        help='output distance/angle restraints')

    parser.add_argument(
        '-t', '--top', required=True,
        action='store', dest='top',
        help='output file for ccg --> cg / fullatom convertation')

    parser.add_argument(
        '-m', '--map', required=True,
        action='store', dest='map',
        help='2d coords map file')

    parser.add_argument(
        '--seq',
        action='store', dest='seq',
        help='scaffold sequence')

    parser.add_argument(
        '--oligs',
        action='store', dest='oligs',
        help='staple sequences exported csv')

    parser.add_argument(
        '--debug',
        action='store_true', dest='debug',
        help='Debug run')

    args = parser.parse_args()

    return args


class COSM(object):

    DEBUG = False

    # Basic lattice properties
    FF = None  # No forcefield by default
    HELIXDIST = None

    # Properties of current project
    NUM = 0
    LEN = 0
    helices = None
    Path = None
    allPath = None
    scaf_linear = True
    scheme = None
    newscheme = None
    insert = None
    delet = None
    mapf = None
    a_sc = 1    # number of particle printed
    nl = 0      # number of PDB chains

    # ------------- output files ----------
    outpdb = []
    outpdbc = []
    outmap = []

    # Unknown arrays. Have to understand their purpose
    scneed = []
    scconn = []
    atomsh = [None]
    anames = [None]
    dif = None
    add_ends = []
    needst = []
    needed = []
    connst = []
    connects = []
    add_tmp = []
    ssoligs = []
    ssind = []
    sstrans1 = []
    sstrans2 = []
    add = []
    PAIR_DONE = []
    R_DONE = []
    T_DONE = []
    pdb = []
    cnct = []

    def __init__(self):
        self.FF = ForceField()

    def parse_json(self, obj):
        # create dictionaries (keyed by virtual helix #) of
        # row/col, scaf array, stap array

        strands = obj['vstrands']

        helices = OrderedDict()

        # Enumerate strands
        for s in strands:
            helices[s['num']] = s

        self.NUM = len(strands)
        self.LEN = len(strands[0]['scaf'])

        # Check if all nums are uniq

        if len(strands) != len(helices.keys()):
            # Get strands number ids
            hNums = [s['num'] for s in strands]
            dups = util.find_duplicates(hNums)
            raise COSMError(
                "USER: "
                "Strans %s appear in input file more than one time" %
                ','.join(map(str, dups))
            )

        self.helices = helices

        trans = []
        trans_sc = []

        scaf_ends = {
            'e5': None,  # 5' end. Should be [helix_num, index],
            'e3': None,  # 3' end. Should be [helix_num, index],
        }

        Path = {'scaf': {}, 'stap': {}}
        allPath = {}

        insert = {}
        deletion = {}

        for num, strand in helices.items():
            strand['dir'] = self.get_strand_direction(strand)

            # - means this cell is empty
            scafPath = ['-'] * self.LEN
            stapPath = ['-'] * self.LEN

            # scaffold path
            for i in range(self.LEN):
                base = strand['scaf'][i]
                # it came from somewhere else
                if (base[0] != num) & (base[0] != -1):
                    scafPath[i] = base[0]
                    trans_sc.append([[num, base[0]]] + [[i, base[1]]])

                # it came from somewhere else
                if (base[2] != num) & (base[2] != -1):
                    scafPath[i] = base[2]
                    trans_sc.append([[num, base[2]]] + [[i, base[3]]])

                # it's part of this helix
                if (base[0] == num) & (base[2] == num):
                    scafPath[i] = num

            # split search to future
            # may be we will have to handle conficts while moving skips
            # if they are on ends
            for i in range(self.LEN):
                base = strand['scaf'][i]

                # Detect 5' end
                if (base[0] == -1) & (base[2] == num):

                    # Check if end can be deleted
                    if strand['skip'][i] != 0:
                        i_ = self.move_skip(strand, i)
                    else:
                        i_ = i

                    scafPath[i] = 'e5'
                    scaf_ends['e5'] = [num, i]

                # Detect 3' end
                if (base[0] == num) & (base[2] == -1):

                    # Check if end can be deleted
                    if strand['skip'][i] != 0:
                        i_ = self.move_skip(strand, i)
                    else:
                        i_ = i

                    scafPath[i] = 'e3'
                    scaf_ends['e3'] = [num, i]

            Path['scaf'][num] = scafPath

            # staple path
            for i in range(self.LEN):
                base = strand['stap'][i]

                # it came from somewhere else
                if (base[0] != num) & (base[0] != -1):
                    stapPath[i] = base[0]
                    trans.append([[num, base[0]]] + [[i, base[1]]])

                # it came from somewhere else
                if (base[2] != num) & (base[2] != -1):
                    stapPath[i] = base[2]
                    trans.append([[num, base[2]]] + [[i, base[3]]])

                # it's part of this helix
                if (base[0] == num) & (base[2] == num):
                    stapPath[i] = num

            # split search to future
            # may be we will have to handle conficts while moving skips
            # if they are on ends
            for i in range(self.LEN):
                base = strand['stap'][i]

                # Detect 5' end
                if (base[0] == -1) & (base[2] == num):
                    # Check if end can be deleted
                    if strand['skip'][i] != 0:
                        i_ = self.move_skip(
                            strand,
                            i,
                            # Because 5' always antiparallel to scaffold
                            -1,
                        )
                    else:
                        i_ = i

                    stapPath[i] = 'e5'

                # Detect 3' end
                if (base[0] == num) & (base[2] == -1):

                    # Check if end can be deleted
                    if strand['skip'][i] != 0:
                        i_ = self.move_skip(strand, i)
                    else:
                        i_ = i

                    stapPath[i] = 'e3'

            Path['stap'][num] = stapPath

            # join scaffold and staple
            joinedPath = ['-'] * self.LEN

            for i in range(self.LEN):
                # if self.DEBUG:
                #     print("Scaf: %s Stap: %s" % (scafPath[i], stapPath[i]))

                if (scafPath[i] == '-') & (stapPath[i] == '-'):
                    # We are already prefilled array with '-'
                    pass
                elif (scafPath[i] != '-') & (stapPath[i] != '-'):
                    # d - means duplex part
                    joinedPath[i] = 'd'
                elif (scafPath[i] != '-') & (stapPath[i] == '-'):
                    # s - means single stranded part
                    joinedPath[i] = 's'
                else:
                    joinedPath[i] = 'o'
                    raise COSMError(
                        "USER: "
                        "Single-stranded staple: helix: %d pos: %d" % (num, i_)
                    )

            allPath[num] = joinedPath

            # ------------- loop --------------

            l_ = np.array(strand['loop'])
            s_ = np.array(strand['skip'])
            # sum loop and skip to remove positions
            # with simultaneous deletion and insertion
            ls_ = l_ + s_
            # keep only positions, where insertions remains
            strand['loop'] = l_ * (ls_ == 1)
            # keep only positions, where deletions remains
            strand['skip'] = s_ * (ls_ == -1)

            # loop corresponds to insertions
            for i in range(self.LEN):
                n = strand['loop'][i]
                # 0 means no deletion or insertion
                if n == 0:
                    pass
                # Normally it should be just 1 nucleotide
                elif n == 1:
                    insert[(num, i)] = n
                # N - means how many nucleotides to insert
                elif n > 1:
                    raise COSMError(
                        'USER: More than one base inserted in one place')
                else:
                    raise COSMError(
                        'USER: Wrong insertion definition')

            # skip corresponds to deletions
            for i in range(self.LEN):
                n = strand['skip'][i]
                # 0 means no deletion or insertion
                if n == 0:
                    pass
                # Normally it should be just 1 nucleotide (-1)
                elif abs(n) == 1:
                    deletion[(num, i)] = abs(n)
                # -N - means how many nucleotides to delete
                elif abs(n) > 1:
                    raise COSMError(
                        'USER: More than one base deleted in one place')
                else:
                    raise COSMError(
                        'USER: Wrong deletion definition')

        if self.check_scaf_ends(scaf_ends):
            # If scaffold is linear ends have to be alredy found
            pass
        else:
            # Ok. Seems that scaffold is circular. Let's check this
            self.scaf_linear = False
            scaf_ends = self.get_circular_scaffold_ends(allPath)
            Path['scaf'][scaf_ends['e5'][0]][scaf_ends['e5'][1]] = 'e5'
            Path['scaf'][scaf_ends['e3'][0]][scaf_ends['e3'][1]] = 'e3'

        self.Path = Path
        self.allPath = allPath
        self.scaf_ends = scaf_ends
        self.trans_sc = trans_sc
        self.trans = trans

    def check_scaf_ends(self, end=None):
        if not end:
            end = self.scaf_ends

        flag = False

        if (
            # Check 5' end
            (isinstance(end['e5'], list) and len(end['e5']) == 2) and
                # Check 3' end
                (isinstance(end['e3'], list) and len(end['e3']) == 2)):

            flag = True

        return flag

    def process(self, args):

        self.DEBUG = args.debug

        # ------------- .json analysis ----------------
        obj = util.load_json(args.input)

        # ------------ lattice ------------
        lattice = util.get_lattice_type(obj)
        self.FF.set_lattice_params(lattice)

        try:
            self.parse_json(obj)
        except COSMError as e:
            raise COSMError(e)
        except:
            if self.DEBUG:
                tb.print_exc()
                raise Exception("Unexpected error:", sys.exc_info())
            else:
                raise COSMError('USER: Unrecoverable error in input json file')

        if self.DEBUG:
            print('Scaf linear:', self.scaf_linear)
            print('ENDS', self.scaf_ends)

    def get_circular_scaffold_ends(self, allPath):
        # ---------------- end5 definition (if circular)  ---------------
        scaf_ends = {
            'e5': None,  # 5' end. Should be [helix_num, index],
            'e3': None,  # 3' end. Should be [helix_num, index],
        }

        dl = 21     # length of duplex self.needed
        if self.FF.SQ:
            dl = 36

        k = min(self.helices.keys())
        end = 0
        while end == 0 and k in allPath:
            # At first look on reversed chains
            # Don't know why though
            if self.helices[k]['dir'] == -1:
                c = self.LEN / 2
                while not end and c < (self.LEN - (dl - 1)):
                    if allPath[k][c:c + dl] == ['d'] * dl:
                        end = c + dl / 2
                    else:
                        c += 1
                c = self.LEN / 2 - 1
                while not end and c >= 0:
                    if allPath[k][c:c + dl] == ['d'] * dl:
                        end = c + dl / 2
                    else:
                        c -= 1
                if not end:
                    k += 1
            else:
                k += 1
            scaf_ends['e5'] = [k, end - 1]
            scaf_ends['e3'] = [k, end]

        # If we didn't find anything
        # look at all chains
        if end == 0:
            k = min(self.helices.keys())
            while end == 0 and k in allPath:
                b = 1
                while not end and b < self.LEN:
                    if allPath[k][b - 1: b + 2] == ['d'] * 3:
                        end = True
                        scaf_ends['e5'] = [k, b]
                        scaf_ends['e3'] = [k, b + self.revers(
                            self.helices[k], 1)]
                    b += 1
                k += 1

        if self.check_scaf_ends(scaf_ends):
            return scaf_ends
        else:
            raise COSMError("ADMIN: PY1. 5'-end search error")

    def load_sequences(self, seqfn, oligfn):

        # ------------- sequence analysis/generation -------------

        seqPath = OrderedDict()
        for num in self.helices.keys():
            seqPath[num] = ['-'] * self.LEN

        if seqfn and oligfn:
            pass

        elif seqfn or oligfn:
            raise COSMError(
                'USER: Both scaffold and staple sequences are self.needed')
        else:
            # Do not know what for this file is
            self.add.append('no seq\n')
            self.seqPath = seqPath
            return False

        # Get sequence

        seq = util.parse_sequence(seqfn)

        oligs = util.parse_staples(oligfn)

        for ol in oligs:

            # This event can not happen due to restrictions of input format
            # if line[2][1] != '?':

            seqPath[ol[0]][ol[1]] = ol[4][0]  # first letter

            # This event can not happen due to restrictions of input format
            # if line[2][1] != '?':

            seqPath[ol[2]][ol[3]] = ol[4][-1]  # last letter

        # Useless check because we've already checked olig
        # for vh in seqPath:
        #     for i in seqPath[vh]:
        #         if i not in compl and i != '-':
        #             raise COSMError('USER: Error in staple csv file')

        seq_ = seq

        if self.scaf_linear:
            self.checkseq(seqPath, seq_)
        else:
            i = 0
            flag = False
            ls = len(seq)

            while i < ls:
                try:
                    self.checkseq(seqPath, seq_)
                    flag = True
                    break
                except COSMError:
                    pass

                i += 1
                seq_ = seq_[1:] + seq_[0]

            if not flag:
                raise COSMError(
                    'USER: Sequence does not match staple list.')

        self.seq = seq_

        # Do not know what for this file is
        self.add.append('seq: ' + self.seq + '\n')

    def checkseq(self, seqPath, seq, start=None, end=None):
        """ Check if sequence and oligs fit on this step """

        if not start:
            start = self.scaf_ends['e5']

        if not end:
            end = self.scaf_ends['e3']

        [vh, base] = start

        i = 0
        ls = len(seq)
        lall = len(self.helices) * self.LEN
        trans = self.trans_sc[:]  # copy array of transitions

        while [vh, base] != self.scaf_ends['e3']:

            if i >= ls:
                raise COSMError(
                    'USER: Sequence to short')
            elif i >= lall:
                raise COSMError(
                    'USER: Cannot reach scaffold end. Code is broken')

            if seqPath[vh][base] != '-':
                o_ = seqPath[vh][base]
                s_ = seq[i]
                if util.compl[s_] != o_:
                    raise COSMError(
                        'USER: Sequence does not match staple list '
                        'at helix: %d position: %d' % (vh, base))

            [vh, base], trans = self.nextbase(vh, base, trans)
            i += 1

    def nextbase(self, vh, base, trans):
        """ Get next base after """
        # works two times in sequence search:
        # 1st iteration for next sequence negin point
        # 2nd iterantion for sequence checking
        k = self.helices[vh]['dir']

        nbase = [vh, base + k]

        for i, pair in enumerate(trans):
            if pair[0][0] == vh and pair[1][0] == base:
                nbase = [pair[0][1], pair[1][1]]

                trans.pop(i)

                # now we need to eliminate bacward transition
                rev = [
                    [pair[0][1], pair[0][0]],
                    [pair[1][1], pair[1][0]]
                ]
                trans.remove(rev)

        return nbase, trans

    def create_segment(self, vh, base):

        seg_ = [self.allPath[vh][base], vh, [base, None]]

        return seg_

    def build_scheme(self, args, seqPath=None, start=None, end=None):
        #
        # ------------ scaffold scheme (path) generation -------------

        if start is None:
            start = self.scaf_ends['e5']

        if end is None:
            end = self.scaf_ends['e3']

        lall = len(self.helices) * self.LEN
        trans = self.trans_sc[:]  # copy array of transitions

        mapf = list()
        scheme = list()

        # start from 5' end
        [vh, base] = start
        seg_ = self.create_segment(vh, base)
        mapf.append('L %d:%d' % (base, vh))

        i = 0

        # scan path
        while [vh, base] != end:

            # backup check for infinite loop
            if i >= lall:
                raise COSMError(
                    'USER: Cannot reach scaffold end. Code is broken')

            # keep previous base
            prev_ = [vh, base]
            # get nextbase
            [vh, base], trans = self.nextbase(vh, base, trans)

            # if we are on the same helix and the same type (duplex or single)
            # go to next base
            if vh == seg_[1] and self.allPath[vh][base] == seg_[0]:
                pass
            else:
                # else - end segment and start new one
                if vh != seg_[1]:
                    mapf.append('L %d:%d' % (prev_[1], prev_[0]))
                    mapf.append('L %d:%d' % (base, vh))

                seg_[2][1] = prev_[1]
                scheme.append(seg_)
                seg_ = self.create_segment(vh, base)

            i += 1

        mapf.append('L %d:%d' % (base, vh))
        seg_[2][1] = base
        scheme.append(seg_)

        self.validate_scheme(scheme)

        self.scheme = scheme
        self.mapf = mapf

    def validate_scheme(self, scheme):
        valPath = copy.deepcopy(self.Path['scaf'])
        valStap = copy.deepcopy(self.Path['stap'])

        # -------- more than one scaffolds validation -------------

        for p in scheme:
            pb = max(p[2])
            pe = min(p[2])
            for i in range(pe, pb + 1):
                valPath[p[1]][i] = '-'

        for vh in valPath:
            for i in valPath[vh]:
                if i != '-':
                    raise COSMError('USER: More than one scaffold chain')
        # ---------------------------------------------------------

        self.valPath = valPath
        self.valStap = valStap

    def build_newscheme(self):
        insert = list()
        deletion = list()

        newscheme = copy.copy(self.scheme)

        for part in newscheme:
            k = self.helices[part[1]]['dir']

            if part[0] == 'd':      # duplexes
                ins1 = []
                ins2 = []
                ins = []
                dele1 = []
                dele2 = []
                dele = []

                beg = part[2][0] - part[2][0] % self.FF.HC + \
                    self.FF.HC * (k + 1) / 2 * bool(part[2][0] % self.FF.HC)
                end = part[2][1] - part[2][1] % self.FF.HC - \
                    self.FF.HC * (k - 1) / 2 * bool(part[2][1] % self.FF.HC)
                isH = False
                for base in range(part[2][0], part[2][1] + k, k):
                    if base % self.FF.HC:
                        isH = True
                if (end - beg) * k >= 0 and isH:
                    part[2] = [(
                        part[2][0], (
                            beg - part[2][0]) * k), range(
                                beg, end + self.FF.HC * k, self.FF.HC * k),
                        (part[2][1], (part[2][1] - end) * k)]

                    # find insertions
                    # first term
                    for i in range(part[2][0][0], part[2][1][0], k):
                        if (part[1], i) in insert:
                            ins1.append(i)

                    part[2][0] = tuple(list(part[2][0]) + [ins1])

                    # last term
                    for i in range(part[2][1][-1] + k, part[2][2][0] + k, k):
                        if (part[1], i) in insert:
                            ins2.append(i)

                    part[2][2] = tuple(list(part[2][2]) + [ins2])

                    if len(part[2][1]) > 1:
                        for i, b in enumerate(part[2][1][:-1]):  # center
                            ins.append([])
                            for l in range(b, b + k * self.FF.HC, k):
                                if (part[1], l) in insert:
                                    ins[-1].append(l)
                        if (part[1], part[2][1][-1]) in insert:
                            ins[-1].append(part[2][1][-1])
                        ins.append([])
                    else:
                        if (part[1], part[2][1][0]) in insert:
                            ins.append([part[2][1][0]])
                        else:
                            ins.append([])
                    part.append(ins)

                    # find deletions

                    for i in range(part[2][0][0], part[2][1][0], k):
                        if (part[1], i) in deletion:
                            dele1.append(i)
                    for i in range(part[2][1][-1] + k, part[2][2][0] + k, k):
                        if (part[1], i) in deletion:
                            dele2.append(i)

                    part[2][0] = tuple(list(part[2][0]) + [dele1])
                    part[2][2] = tuple(list(part[2][2]) + [dele2])

                    if len(part[2][1]) > 1:
                        for i, b in enumerate(part[2][1][:-1]):
                            dele.append([])
                            for l in range(b, b + k * self.FF.HC, k):
                                if (part[1], l) in deletion:
                                    dele[-1].append(l)
                        if (part[1], part[2][1][-1]) in deletion:
                            dele[-1].append(part[2][1][-1])
                        dele.append([])
                    else:
                        if (part[1], part[2][1][0]) in deletion:
                            dele.append([part[2][1][0]])
                        else:
                            dele.append([])
                    part.append(dele)

                    if not part[2][0][1]:
                        part[2][0] = None
                    if not part[2][2][1]:
                        part[2][2] = None

                else:      # if duplex is small and without H - N particles
                    part[2] = range(part[2][0], part[2][1] + k, k)
                    part[0] = 'n'
                    for i in range(part[2][0], part[2][-1] + k, k):
                        if (part[1], i) in insert:
                            ins.append(i)
                        if (part[1], i) in deletion:
                            dele.append(i)
                    part.append(ins)
                    part.append(dele)
            elif part[0] == 's':    # single-stranded
                part[2] = range(part[2][0], part[2][1] + k, k)
                for i in part[2]:
                    if (part[1], i) in insert or (part[1], i) in deletion:
                        raise COSMError(
                            'USER: Insertion or deletion '
                            'in single-stranded region'
                        )

        self.insert = insert
        self.deletion = deletion
        self.newscheme = newscheme

    def assign_ff_particles(self):
        for a, part in enumerate(self.newscheme):
            if part[0] == 'd':
                # first term

                if part[2][0]:
                    n = 'T' + str(part[2][0][1])
                m = 'T'
                tmpins = 0
                tmpdel = 0
                if part[2][0]:
                    self.createatom(
                        part[2][0][0],
                        part[1], n,
                        part[2][0][1], 'b',
                        ('t', part[2][0][2], part[2][0][3])
                    )

                    if len(part[2][1]) <= 1:
                        if part[2][2]:
                            self.createatom(
                                part[2][1][0],
                                part[1], 'PT',
                                1, 0,
                                (part[3][0], part[4][0])
                            )

                        else:
                            self.createatom(
                                part[2][1][0],
                                part[1], 'T',
                                1, 0,
                                (part[3][0], part[4][0])
                            )

                    else:
                        self.createatom(
                            part[2][1][0],
                            part[1], 'PT',
                            self.FF.HC, 0,
                            (part[3][0], part[4][0])
                        )

                else:
                    self.createatom(
                        part[2][1][0],
                        part[1], m,
                        self.FF.HC, 'b',
                        (part[3][0], part[4][0])
                    )

                # center
                for i, atom in enumerate(part[2][1][1:-1]):
                    self.createatom(
                        atom,
                        part[1], 'H',
                        self.FF.HC, 0,
                        (part[3][i + 1], part[4][i + 1])
                    )

                # last term
                if part[2][2]:
                    if len(part[2][1]) > 1:
                        self.createatom(
                            part[2][1][-1],
                            part[1], 'PT',
                            part[2][2][1], 0,
                            (part[3][-1], part[4][-1])
                        )

                    self.createatom(
                        part[2][2][0],
                        part[1], 'T' + str(part[2][2][1]),
                        1, 'e',
                        ('e', part[2][2][2], part[2][2][3])
                    )

                else:
                    if len(part[2][1]) > 1:
                        self.createatom(
                            part[2][1][-1],
                            part[1], 'T',
                            1, 'e',
                            (part[3][-1], part[4][-1])
                        )

            elif part[0] == 's':
                if a == 0:
                    self.createatom(
                        part[2][0],
                        part[1], 'S',
                        1, 's',
                        False
                    )

                    if len(part[2]) > 1:
                        for atom in part[2][1:]:
                            self.createatom(
                                atom,
                                part[1], 'S',
                                1, 's',
                                False
                            )

                else:
                    for atom in part[2]:
                        self.createatom(
                            atom,
                            part[1], 'S',
                            1, 's',
                            False
                        )

            elif part[0] == 'n':
                tmpins = False
                tmpdel = False
                if len(part[2]) <= 2:
                    for atom in part[2]:
                        if atom in part[4]:
                            raise COSMError(
                                'USER: Deletion at terminal base pair')
                        if atom in part[3]:
                            self.createatom(
                                atom,
                                part[1], 'T',
                                1, 0,
                                (1, 0)
                            )

                        else:
                            self.createatom(
                                atom,
                                part[1], 'T',
                                1, 0,
                                (0, 0)
                            )

                else:

                    if part[2][0] in part[4]:
                        raise COSMError('USER: Deletion at terminal particle')

                    if part[2][-1] in part[4]:
                        raise COSMError('USER: Deletion at terminal particle')

                    if part[2][0] in part[3]:
                        self.createatom(
                            part[2][0],
                            part[1], 'T',
                            1, 0,
                            (1, 0)
                        )

                    else:
                        self.createatom(
                            part[2][0],
                            part[1], 'T',
                            1, 0,
                            (0, 0)
                        )

                    for atom in part[2][1:-1]:
                        tmpins = 0
                        tmpdel = 0
                        if atom in part[3]:
                            tmpins = 1
                        if atom in part[4]:
                            tmpdel = 1
                        self.createatom(
                            atom,
                            part[1], 'N',
                            1, 0,
                            (tmpins, tmpdel)
                        )

                    if part[2][-1] in part[3]:
                        self.createatom(
                            part[2][-1],
                            part[1], 'T',
                            1, 0,
                            (1, 0)
                        )

                    else:
                        self.createatom(
                            part[2][-1],
                            part[1], 'T',
                            1, 0,
                            (0, 0)
                        )

            else:
                raise COSMError("ADMIN: PY1. Unknown chain type")

    def createatom(self, z, chain, tname, length, end, modif):
        """Print particles in the interval -- calculate coords etc"""
        nls = string.ascii_uppercase + string.digits    # for additional chains
        nlname = nls[self.nl]

        row = self.helices[chain]['row']
        col = self.helices[chain]['col']

        if not self.FF.SQ:
            ox, oy = self.honeycomb(row, col)
        else:
            ox, oy = self.square(row, col)
        if self.a_sc == 1:
            self.dif = None
        aname = 'D'
        if modif:
            if modif[0] not in ['t', 'e']:
                inserts = modif[0]
                delet = modif[1]
                instype = None
            else:
                instype = modif[0]
                inserts = modif[1]
                delet = modif[2]
        else:
            inserts = None
            delet = None
        if inserts or delet:
            if isinstance(inserts, int):
                inserts = [inserts]
            if isinstance(delet, int):
                delet = [delet]
            if set(inserts) & set(delet):
                raise COSMError(
                    'USER: Insertion and deletion in the same place.')
        k = self.helices[chain]['dir']

# THIS PART CAN NEVER BEEN EXECUTED
# BECAUSE self.unknown_func_1 can never been finised
# due to exception in the very beggining
#        if (chain, z) in self.ssnear:
#            n = self.ssnear.index((chain, z))
#            self.ssnear[n] = self.a_sc
#        elif (chain, z - 1) in self.ssnear:
#            n = self.ssnear.index((chain, z - 1))
#            self.ssnear[n] = self.a_sc
#        if (chain, z) in self.ssoligs:
#            n = self.ssoligs.index((chain, z))
#            if self.ssnear[n]:
#                self.add.append(
#                    'scaf ; ' + str(self.ssnear[n]) +
#                    ' ; ' + str(self.a_sc) + '\n'
#                    )
        savedif = False
        if tname == 'N':
            if delet:
                self.addtoout(aname, 'N', nlname, ox, oy, z, chain, 1, 'D')
            if inserts:
                self.addtoout(aname, 'N', nlname, ox, oy, z, chain, 1, 'I')
        else:
            if self.dif and tname not in self.FF.term:
                tname = 'PT'
            if self.dif and tname in self.FF.term and tname[0] == 'T':
                if not inserts and not delet:
                    tname = 'T'
                else:
                    length = int(tname[1])
                    tname = 'T'
                    savedif = True
            self.dif = self.findconnects(chain, z, self.a_sc, length, end)
            if self.dif:
                if tname not in self.FF.term:
                    tname = 'PT'
                if tname in self.FF.term and tname[0] == 'T':
                    tname = 'T'
                if z in inserts:
                    M = 'I'
                elif z in delet:
                    M = 'D'
                else:
                    M = None
                self.addtoout(aname, tname, nlname, ox, oy, z, chain, 1, M)
                for i in range(1, length):
                    if z + i * k in inserts:
                        M = 'I'
                    elif z + i * k in delet:
                        M = 'D'
                    else:
                        M = None
                    self.addtoout(
                        aname, 'N', nlname, ox, oy, z + i * k, chain, 1, M)
            elif inserts or delet:
                if not instype:
                    if tname not in self.FF.term:
                        tname = 'PT'
                    if tname in self.FF.term and tname[0] == 'T':
                        tname = 'T'
                    if z in inserts:
                        M = 'I'
                    elif z in delet:
                        M = 'D'
                    else:
                        M = None
                    if length > 1:
                        self.addtoout(
                            aname, tname, nlname, ox, oy, z, chain, length, M)
                    elif M:
                        self.addtoout(
                            aname, tname, nlname, ox, oy,
                            z, chain, length, M + 'e')  # first insert
                    else:
                        self.addtoout(
                            aname, tname, nlname,
                            ox, oy, z,
                            chain, length, None
                        )

                    for i in range(1, length):
                        if z + i * k in inserts:
                            M = 'I'
                        elif z + i * k in delet:
                            M = 'D'
                        else:
                            M = None
                        self.addtoout(
                            aname, 'N', nlname, ox, oy, z + i * k, chain, 1, M)
                    self.dif = True
                else:
                    if not savedif:
                        length = self.FF.count[tname]
                    if z in delet and tname == 'T1':
                        # because of T1
                        raise COSMError(
                            'USER: Deletion at terminal base pair')
                        MT = 'D'
                    if z in inserts:
                        MT = 'I'
                    else:
                        MT = None
                    if instype == 't':
                        self.addtoout(
                            aname, 'T', nlname, ox, oy, z, chain, length, MT)
                        for i in range(1, length):
                            M = None
                            if z + i * k in inserts:
                                M = 'I'
                            elif z + i * k in delet:
                                M = 'D'
                            self.addtoout(
                                aname, 'N', nlname,
                                ox, oy, z + i * k,
                                chain, 1, M
                            )
                    else:
                        for i in range(1, length):
                            M = None
                            if z + (i - length) * k in inserts:
                                M = 'I'
                            elif z + (i - length) * k in delet:
                                M = 'D'
                            self.addtoout(aname, 'N', nlname, ox, oy, z + (
                                i - length) * k, chain, 1, M)
                        if MT:
                            self.addtoout(
                                aname, 'T', nlname, ox, oy, z,
                                chain, length, MT + 'e')
                        else:
                            self.addtoout(
                                aname, 'T', nlname, ox, oy, z,
                                chain, 1, None)  # !!!
            else:
                self.addtoout(
                    aname, tname, nlname, ox, oy, z,
                    chain, length, None)  # quite normal atom

    def findconnects(self, vh, j, number, length, end):
        """Find crossover restraints"""
        if end == 's':
            return None
        # [x, y] = [0, 0]
        k = self.helices[vh]['dir']
        seq = [j]
        seq += range(j + k, j + k * length, k)
        ind = []
        for i in seq:
            dif = (i - j) * k
            if self.Path['stap'][vh][i] not in [vh, 'e5', 'e3', '-']:
                if [vh, i] in self.needst:
                    n = self.needst.index([vh, i])
                    if i == j:
                        self.needst[n] = (number, 0)
                        self.add_tmp.append(
                            [number, 0, self.connst[n][0], self.connst[n][1]])
                    elif i % self.FF.HC == self.FF.HC - 1:
                        self.needst[n] = (number, dif)
                        if k == -1:
                            self.add_tmp.append(
                                [
                                    number, 1, self.connst[n][0],
                                    self.connst[n][1]
                                ])
                        else:
                            self.add_tmp.append(
                                [
                                    number,
                                    length - 1,
                                    self.connst[n][0],
                                    self.connst[n][1]
                                ]
                            )
                    else:
                        self.needst[n] = (number + dif, 0)
                        self.add_tmp.append(
                            [
                                number + dif,
                                0,
                                self.connst[n][0],
                                self.connst[n][1]
                            ])

                if [vh, i] in self.needed:
                    n = self.needed.index([vh, i])
                    if i == j:
                        self.needed[n] = number
                    elif i % self.FF.HC == self.FF.HC - 1:
                        if k == 1:
                            self.needed[n] = number + 1
                        else:
                            self.needed[n] = number
                    else:
                        ind.append(dif)
                        self.needed[n] = number + dif
        for i in seq:
            dif = (i - j) * k
            if self.Path['stap'][vh][i] not in [vh, 'e5', 'e3', '-']:
                if [vh, i] not in self.needed:
                    l = self.Path['stap'][vh][i]
                    for pair in self.trans:
                        if pair[0] == [vh, l] and pair[1][0] == i:
                            self.needst.append([pair[0][1], pair[1][1]])
                            if i == j:
                                self.connst.append((number, 0))
                                if number not in self.connects + self.needed:
                                    self.connects.append(number)
                                    self.needed.append(
                                        [pair[0][1], pair[1][1]])
                            elif i % self.FF.HC == self.FF.HC - 1:
                                if k == -1:
                                    self.connst.append((number, 1))
                                else:
                                    self.connst.append((number, length - 1))
                                if [vh, i + 1] not in self.needed:
                                    if k == 1:
                                        if (
                                            number + 1 not in
                                                self.connects + self.needed):

                                            self.connects.append(number + 1)
                                            self.needed.append(
                                                [pair[0][1], pair[1][1]])
                                    else:
                                        if (
                                            number not in
                                                self.connects + self.needed):

                                            self.connects.append(number)
                                            self.needed.append(
                                                [pair[0][1], pair[1][1]])
                            else:
                                self.connst.append((number + dif, 0))
                                if self.dif not in ind:
                                    ind.append(dif)
                                self.connects.append(number + dif)
                                self.needed.append([pair[0][1], pair[1][1]])

        if ind:
            return ind
        return None

    def process_next_(self):

        self.cnct.append('[ distance_restraints ]\n')
        self.cnct.append('; staple crossovers\n')

        l = 1

        # -------- check circular staples --------

        for vh in self.valStap:
            s_ = self.valStap[vh]
            for i_ in range(self.LEN):
                if s_[i_] != '-':
                    raise COSMError(
                        'USER: Circular staple: helix: %d pos: %d' % (vh, i_))

        # ---------------------------------------

        done = []
        num = len(self.outpdb) + 1
        while len(done) < len(self.ssoligs):
            self.outpdb.append(['TER\n'])
            self.nl += 1
            i = 0
            while self.ssoligs[i] in done or self.ssind[i] != 'b':
                i += 1
            now = self.ssoligs[i]
            done.append(now)
            [vh, base] = self.ssoligs[i]
            k = self.helices[vh]['dir']
            self.createatom(base, vh, 'OT', 1, 's')
            n = self.ssoligs.index((vh, base))
            if self.ssnear[n]:
                self.cnct.append(
                    str(self.ssnear[n]) + '\t' +
                    str(self.a_sc - 1) + '\t' +
                    '1\t' +
                    str(l) +
                    '\t1\t0.33\t0.35\t0.37\t0.2\t; ssnear\n')

                l += 1
                self.add.append(
                    'cross ; ' +
                    ' '.join(map(str, [self.ssnear[n], 0])) +
                    ' ' +
                    ' '.join(map(str, [self.a_sc - 1, 0])) +
                    '\n'
                    )

            base += k
            while (vh, base) not in self.ssoligs:
                k = self.helices[vh]['dir']
                self.createatom(base, vh, 'O', 1, 's')
                if (vh, base) in self.sstrans1 and (vh, base) not in done:
                    n = self.sstrans1.index((vh, base))
                    (vh, base) = self.sstrans2[n]
                    done.append((vh, base))
                elif (vh, base) in self.sstrans2 and (vh, base) not in done:
                    n = self.sstrans2.index((vh, base))
                    (vh, base) = self.sstrans1[n]
                    done.append((vh, base))
                else:
                    base += k
            done.append((vh, base))
            self.createatom(base, vh, 'O', 1, 's')
            n = self.ssoligs.index((vh, base))
            if self.ssnear[n]:
                self.cnct.append(
                    str(self.ssnear[n]) + '\t' +
                    str(self.a_sc - 1) + '\t' +
                    '1\t' +
                    str(l) +
                    '\t1\t0.33\t0.35\t0.37\t0.2\n'
                    )

                l += 1
                self.add.append(
                    'cross ; ' +
                    ' '.join(map(str, [self.ssnear[n], 0])) +
                    ' ; ' +
                    ' '.join(map(str, [self.a_sc - 1, 0])) +
                    '\n'
                    )

        self.outpdb[0][3] += 'T'

        for i, c in enumerate(self.atomsh):
            if i:
                self.allPath[c[0]][c[1]] = i
        for vh in self.allPath:
            k = self.helices[vh]['dir']
            if k == -1:
                self.allPath[vh] = self.allPath[vh][::-1]
            for i, c in enumerate(self.allPath[vh]):
                if c not in ['-', 'd', 's', 'DELE']:
                    ntmp = c
                    ctmp = 0
                    self.allPath[vh][i] = (c, 0)
                elif c in ['d', 's']:
                    ctmp += 1
                    self.allPath[vh][i] = (ntmp, ctmp)
            if k == -1:
                self.allPath[vh] = self.allPath[vh][::-1]

        for i in self.add_ends:
            p_i = self.allPath[i[1]][i[2]]
            if p_i == 'DELE':
                p_i = self.allPath[i[1]][i[2] - self.helices[vh]['dir']]
            self.add.append(
                i[0] + ' ; ' +
                str(p_i[0]) + ' ' +
                str(p_i[1]) +
                '\n'
                )

        for pair in self.trans:
            [[b1, b2], [a1, a2]] = pair
            if not [[b2, b1], [a2, a1]] in self.PAIR_DONE:
                self.PAIR_DONE.append(pair)
                if not (a1 % self.FF.HC or a2 % self.FF.HC):
                    i1 = self.nsearch(b1, a1, 1)
                    i2 = self.nsearch(b2, a2, 1)
                    self.wr_add(i1, i2)
                    self.wr_restr(i1, i2, l)
                    if [[b1, b2], [a1 - 1, a2 - 1]] in self.trans:
                        i1 = self.nsearch(b1, a1 - 1, -1)
                        i2 = self.nsearch(b2, a2 - 1, -1)
                        self.wr_add(i1, i2)
                elif (
                    a1 % self.FF.HC == self.FF.HC - 1 and
                        a2 % self.FF.HC == self.FF.HC - 1):

                    if [[b1, b2], [a1 + 1, a2 + 1]] not in self.trans:
                        i1 = self.nsearch(b1, a1, -1)
                        i2 = self.nsearch(b2, a2, -1)
                        self.wr_add(i1, i2)
                        self.wr_restr(i1, i2, l)
                else:
                    i1 = self.nsearch(b1, a1, 1)
                    i2 = self.nsearch(b2, a2, 1)
                    self.wr_add(i1, i2)
                    self.wr_restr(i1, i2, l)

        prevB = False
        for atom in self.outpdb:
            if len(atom) > 1:
                if (atom[3][0] == 'T' and
                    atom[1] < self.a_sc - 1 and
                    atom[1] > 1 and
                    len(self.outpdb[atom[1]]) > 1 and
                        not prevB):

                    if (
                        self.outpdb[atom[1]][3] != 'S' and
                            self.outpdb[atom[1] - 2][3] != 'S'):

                        new_name = self.TtoB(atom[3], atom[1])
                        prevB = True
                        self.outpdb[atom[1] - 1][3] = new_name
                        self.outmap[atom[1] - 1][1] = new_name
                        if new_name[0] == 'B':
                            if self.outpdb[atom[1]][3] == 'T':
                                self.outpdb[atom[1]][3] = 'B'
                                self.outmap[atom[1]][1] = 'B'
                            elif self.outpdb[atom[1]][3][0] == 'T':
                                self.outpdb[atom[1]][
                                    3] = 'B' + self.outpdb[atom[1]][3][1:]
                                self.outmap[atom[1]][
                                    1] = 'B' + self.outmap[atom[1]][1][1:]
                elif prevB:
                    prevB = False
            self.outatom(atom)
        for atom in self.outmap:
            # atom[3] = vhNums.index(atom[3])
            self.mapf.append(
                "{0[0]:s} {0[1]:s}:{0[2]:d}:{0[3]:d}:{0[4]:d}\n".format(atom))

        self.cnct.append('; scaffold crossovers\n')
        if not self.scaf_linear:
            self.cnct.append(
                '1\t' +
                str(num - 1) +
                '\t1\t' +
                str(l) +
                '\t1\t0.33\t0.35\t0.37\t2.5\t; 5\' - 3\' ends' +
                '\n')
        l += 1

        for i in range(len(self.scconn)):
            if not (
                isinstance(self.scconn[i], int) and
                isinstance(self.scneed[i], int)
                    ):

                raise COSMError('ADMIN: PY1. Error with scaffold restraints')

            self.cnct.append(
                str(self.scconn[i]) + '\t' +
                str(self.scneed[i]) + '\t' +
                '1\t' +
                str(l) +
                '\t1\t0.33\t0.35\t0.37\t2.5\n'
                )

            l += 1
            self.add.append(
                'scaf ; ' +
                str(self.scconn[i]) +
                ' ; ' + str(self.scneed[i]) +
                '\n'
                )
        self.cnct.append('\n')

        self.cnct.append('; lattice restraints\n')
        if len(self.connects) != len(self.needed):
            print('Smth with self.connects is wrong')  # ~ error
        self.cnct.append('\n')

        template = "{0[0]:<6s}{0[1]:>5d}{0[2]:>5d}"
        for i in range(1, num):     # doesn't work for ss oligs
            self.outpdbc.append(template.format(['CONECT', i, i + 1]) + '\n')

        for conect in self.outpdbc:
            self.pdb.append(conect)

        # ---------- lattice restraints ---------

        name = 'H'

        if self.FF.SQ:
            [a, b, c] = ['4.34', '4.4', '4.43']
        else:
            [a, b, c] = ['3.64', '3.81', '4.16']

        triples = {}
        trcoords = {}
        for i, atom in enumerate(self.atomsh):
            if atom and self.anames[i] in [name, 'T', 'TT', 'PT']:
                x = atom[1]
                trcoords[x] = trcoords.get(x, []) + [
                    (
                        self.helices[atom[0]]['row'],
                        self.helices[atom[0]]['col']
                    )]
                triples[x] = triples.get(x, []) + [i]
        for x in trcoords:
            vhs = trcoords[x]
            for i, base in enumerate(vhs):
                tr = self.findtriples(base, vhs)
                for pair in tr:
                    one = self.atomsh[triples[x][i]]
                    two = self.atomsh[triples[x][pair[0]]]
                    thr = self.atomsh[triples[x][pair[1]]]
                    if self.crosscheck(one, two, thr):
                        self.cnct.append(
                            str(triples[x][i]) + '\t' +
                            str(triples[x][pair[1]]) + '\t' +
                            '1\t' +
                            str(l) + '\t' +
                            '1\t' + a + '\t' + b + '\t' + c +
                            '\t1.5' +
                            '\n')

                        l += 1

    @staticmethod
    def get_strand_direction(s):
        """Determine direction of scaffold or staple chain (1/-1)"""
        if (s['row'] + s['col']) % 2 == 0:
            d = 1
        else:
            d = -1

        return d

    def checkreg(self, pair):
        [[vh1, b1], [vh2, b2]] = pair
        if not (isinstance(b1, int) and isinstance(b2, int)):
            return True
        c = -1
        ok1 = False
        ok2 = False
        end = ['-']
        if self.scaf_linear:
            end += ['e5', 'e3']
        while (self.Path['scaf'][vh1][b1 + c] not in end and
                self.Path['scaf'][vh2][b2 + c] not in end and
                c > self.FF.LIM * -1):
            if self.Path['stap'][vh1][b1 + c] == vh2:
                ok1 = True
            c -= 1
        c = 1

        while (self.Path['scaf'][vh1][b1 + c] not in end and
                self.Path['scaf'][vh2][b2 + c] not in end and
                c < self.FF.LIM):
            if self.Path['stap'][vh1][b1 + c] == vh2:
                ok2 = True
            c += 1
        return ok1 | ok2

    def crosscheck(self, one, two, thr):
        '''Check if crossovers in the range'''
        ok = True
        for pair in [[one, two], [two, thr]]:
            ok = ok & self.checkreg(pair)
        return ok

    def TtoB(self, name, atom1):
        if atom1 == self.a_sc - 1 or atom1 == 1:
            return name
        atom2 = atom1 + 1
        atom1 = self.atomsh[atom1]
        atom2 = self.atomsh[atom2]
        if self.checkreg([atom1, atom2]):
            return name
        else:
            if name == 'T':
                return 'B'
            else:
                return 'B' + name[1:]

    def check_circular(self, vh, i):
        end5 = (vh, i)
        end3 = False
        nextbase = None
        prevcr = False
        while not end3 and nextbase != end5:
            self.valStap[vh][i] = '-'
            if nextbase:
                (vh, i) = nextbase
            else:
                (vh, i) = end5
            k = self.helices[vh]['dir'] * -1  # staples
            nextbase = None
            if self.Path['stap'][vh][i] in [vh, 'e5'] or prevcr:
                nextbase = (vh, i + k)
                prevcr = False
            else:
                if not prevcr:
                    for pair in self.trans:
                        if pair[0][0] == vh and pair[1][0] == i:
                            nextbase = (pair[0][1], pair[1][1])
                            prevcr = True
                if not nextbase:
                    if self.Path['stap'][vh][i] == 'e3':
                        end3 = True
                        self.valStap[vh][i] = '-'
                    else:
                        raise COSMError(
                            'ADMIN: Smth wrong with staple crossover paths')
        if nextbase == end5:
            raise COSMError(
                'ADMIN: PY1. Smth wrong with staple crossover paths')

    def staplesends(self, vh, i, number, length):
        """Find staple ends through path"""
        k = self.helices[vh]['dir']
        if not isinstance(i, int):
            return False
        l = 0
        while l < length:
            a = self.Path['stap'][vh][i]
            if a == 'e5':
                self.add_ends.append(('e5', vh, i))
                self.check_circular(vh, i)
            elif a == 'e3':
                self.add_ends.append(('e3', vh, i))
            i += k
            l += 1

    def scaffoldcross(self, vh, i, number):
        """Find and print scaffold crossovers"""
        k = self.helices[vh]['dir']

        if (vh, i) in self.scneed:  # test this
            n = self.scneed.index((vh, i))
    #        t = self.atomsh[self.scconn][n]
    #        if number
            self.scneed[n] = number
        elif (vh, i) in self.scconn:
            n = self.scconn.index((vh, i))
            self.scconn[n] = number
        else:
            n = None
            for pair in self.trans_sc:
                if pair[0][0] == vh and pair[1][0] == i:
                    if self.Path['scaf'][vh][i + k] not in [vh, '-']:
                        if self.Path['scaf'][vh][i - k] in [vh, '-']:
                            if (vh, i + k) in self.atomsh:
                                n = self.atomsh.index((vh, i + k))
                            if n != number - 1:
                                self.scneed.append((vh, i + k))
                                self.scconn.append(number)
                        else:
                            raise COSMError(
                                'ADMIN: PY1. Error with scaffold '
                                'restraints (2bp duplex)')
                    elif self.Path['scaf'][vh][i - k] not in [vh, '-']:
                        if (vh, i - k) in self.atomsh:
                            n = self.atomsh.index((vh, i - k))
                        if n != number - 1:
                            self.scconn.append((vh, i - k))
                            self.scneed.append(number)

    def outatom(self, atom):
        """Print atom in pdb"""
        if len(atom) > 1:
            template = ("{0[0]:<6s}{0[1]:>5d}  {0[2]:<4s}{0[3]:<3s} "
                        "{0[4]:>1s}{0[5]:>4d}    "
                        "{0[6]:>8.3f}{0[7]:>8.3f}{0[8]:>8.3f}"
                        "{0[9]:>6.2f}{0[10]:>6.2f}            {0[11]:<2s}"
                        )
            t = tuple(atom[:-2])
            self.pdb.append(template.format(t) + '\n')
        else:
            self.pdb.append(atom[0])

    def addtoout(self, aname, tname, nlname, ox, oy, z, chain, length, mod):
        if tname:
            t1 = [aname] + [tname] + [nlname]
            t2 = [ox, oy] + [z * 3.4] + [1.00, 0.00] + [aname, chain, z]
            if mod == 'I':
                self.outpdb.append(['ATOM', self.a_sc] + t1 + [self.a_sc] + t2)
                self.outmap.append(['T', tname, self.a_sc, chain, z])
                self.staplesends(chain, z, self.a_sc, length)
                self.scaffoldcross(chain, z, self.a_sc)
                self.atomsh.append((chain, z))
                self.anames.append(tname)
                self.a_sc += 1
                t1[1] = 'N'
                self.outpdb.append(['ATOM', self.a_sc] + t1 + [self.a_sc] + t2)
                self.outmap.append(['I', 'N', self.a_sc, chain, z])

                self.atomsh.append((chain, z))
                self.anames.append('N')
                self.a_sc += 1
            elif mod == 'Ie':
                savename = tname
                self.outpdb.append(['ATOM', self.a_sc] + t1 + [self.a_sc] + t2)
                self.outmap.append(['I', tname, self.a_sc, chain, z])
                self.atomsh.append((chain, z))
                self.anames.append('N')
                self.a_sc += 1
                t1[1] = savename
                self.outpdb.append(['ATOM', self.a_sc] + t1 + [self.a_sc] + t2)
                self.outmap.append(['T', tname, self.a_sc, chain, z])
                self.staplesends(chain, z, self.a_sc, length)
                self.scaffoldcross(chain, z, self.a_sc)
                self.atomsh.append((chain, z))
                self.anames.append(savename)
                self.a_sc += 1
            elif not mod:
                self.outpdb.append(['ATOM', self.a_sc] + t1 + [self.a_sc] + t2)
                self.outmap.append(['T', tname, self.a_sc, chain, z])
                self.staplesends(chain, z, self.a_sc, length)
                self.scaffoldcross(chain, z, self.a_sc)
                self.atomsh.append((chain, z))
                self.anames.append(tname)
                self.a_sc += 1
            elif mod == 'D':
                self.allPath[chain][z] = 'DELE'
            else:
                raise COSMError('ADMIN: PY1. Wrong modification type')

    def findtriples(self, base, coords):
        (row, col) = base
        res = []
        cpair = [None, None]
        rpair1 = [None, None]
        rpair2 = [None, None]
        rpair3 = [None, None]
        rpair4 = [None, None]
        if self.FF.SQ:
            for i, atom in enumerate(coords):
                if atom == (row, col + 1):
                    cpair[0] = i
                elif atom == (row, col + 2):
                    cpair[1] = i
                elif atom == (row + 1, col):
                    rpair1[0] = i
                elif atom == (row + 2, col):
                    rpair1[1] = i
        else:
            for i, atom in enumerate(coords):
                if atom == (row, col + 1):
                    cpair[0] = i
                elif atom == (row, col + 2):
                    cpair[1] = i
                if (row + col) % 2:
                    if atom == (row + 1, col):
                        rpair3[0] = i
                        rpair4[0] = i
                    elif atom == (row + 1, col + 1):
                        rpair3[1] = i
                    elif atom == (row + 1, col - 1):
                        rpair4[1] = i
                else:
                    if atom == (row, col + 1):
                        rpair1[0] = i
                    elif atom == (row + 1, col + 1):
                        rpair1[1] = i
                    elif atom == (row, col - 1):
                        rpair2[0] = i
                    elif atom == (row + 1, col - 1):
                        rpair2[1] = i
        for pair in [cpair, rpair1, rpair2, rpair3, rpair4]:
            if None not in pair:
                res += [pair]
        return res

    def move_skip(self,
                  strand,
                  i,
                  rev=1  # 1 - direct, -1 - reverse
                  ):
            # Potential problem.
            # May be check boundaries
        skip = strand['skip']
        d = strand['dir'] * rev

        i_ = i
        while skip[i_] != 0:
            i_ += d

        # zero moved cells
        if i_ > i:
            for t in range(i_, i, -1):
                skip[t] = skip[t - 1]
            skip[i] = 0

        else:
            for t in range(i_, i,):
                skip[t] = skip[t + 1]
            skip[i] = 0

        return i_ - i

    def nsearch(self, b, a, direct):
        # determine direct
        i = self.allPath[b][a]
        while i == 'DELE':
            i = self.nsearch(b, a + direct, direct)
        return i

    def wr_restr(self, i1, i2, l):
        if (
            i1[1] not in [0, 1, self.FF.HC - 1] or
                i2[1] not in [0, 1, self.FF.HC - 1]):
            raise COSMError('ADMIN: PY1. Restraints error')

        if (
            (i1[0], i2[0]) not in self.R_DONE and
                (i2[0], i1[0]) not in self.R_DONE):
            self.cnct.append(
                str(i1[0]) + '\t' +
                str(i2[0]) + '\t' +
                '1\t' +
                str(l) +
                '\t1\t1.8\t1.85\t1.9\t1.4' +
                '\n')
            l += 1
            self.R_DONE.append((i1[0], i2[0]))
            self.prconnect(i1[0], i2[0])

    def wr_add(self, i1, i2):
        if (i1, i2) not in self.T_DONE and (i1, i2) not in self.T_DONE:
            self.add.append(
                'cross ; ' + ' '.join(map(str, i1)) + ' ' +
                ' '.join(map(str, i2)) + '\n')
            self.T_DONE.append((i1, i2))

    def lat2xy(self):
        pass

    def honeycomb(self, row, col, dist=None):
        """Cadnano honeycomb lattice into cartesian"""
        # k = distance between helixes
        if dist is None:
            dist = self.FF.HELIXDIST
        ANGLE = 30.0
        x = col * dist * math.cos(math.radians(ANGLE))
        if (row % 2) != 0:
            y = (row * 3 + 1 - col % 2) * \
                dist * math.sin(math.radians(ANGLE))
        else:
            y = (row * 3 + col % 2) * dist * math.sin(math.radians(ANGLE))
        return x, y

    def square(self, row, col, dist=None):
        """Cadnano square lattice into cartesian"""

        if dist is None:
            dist = self.FF.HELIXDIST

        x = col * dist
        y = row * dist
        return x, y

    def prconnect(self, base1, base2):
        """Print crossover restraints"""
        template = "{0[0]:<6s}{0[1]:>5d}{0[2]:>5d}"
        t = tuple(['CONECT', base1, base2])
        self.outpdbc.append(template.format(t) + '\n')
        if not (isinstance(base2, int) and isinstance(base1, int)):
            raise COSMError('ADMIN: PY1. Error with staple restraints')

    def unknown_func_1(self):
        # ------------ print staple crossovers ------------
        for vh in self.allPath.keys():
            h = self.allPath[vh]
            k = self.helices[vh]['dir']
            d = 0
            for i, base in enumerate(h):
                if base == 'o':
                    raise COSMError('USER: Single-stranded staple')
                    if self.Path['stap'][vh][i] in ['e3', 'e5']:
                        self.ssoligs.append((vh, i))
                        self.ssnear.append(None)
                        if (d and k == 1) or (not d and k == -1):
                            self.ssind.append('e')
                        else:
                            self.ssind.append('b')
                        d = not d
                    elif self.Path['stap'][vh][i] == vh:
                        w = False
                        if self.allPath[vh][i - k] == 'd':
                            self.ssoligs.append((vh, i))
                            self.ssnear.append((vh, i - k))
                            w = True
                        elif self.allPath[vh][i + k] == 'd':
                            self.ssoligs.append((vh, i))
                            self.ssnear.append((vh, i + k))
                            w = True
                        if w:
                            if (d and k == 1) or (not d and k == -1):
                                self.ssind.append('e')
                            else:
                                self.ssind.append('b')
                            d = not d
                    else:
                        v1 = self.Path['stap'][vh][i]
                        for pair in self.trans:
                            if (
                                pair[0] in [[vh, v1], [v1, vh]] and
                                    i == pair[1][0]):

                                i1 = pair[1][1]
                                if self.allPath[v1][i1] == 'd':
                                    self.ssoligs.append((vh, i))
                                    self.ssnear.append((v1, i1))
                                    if (d and k == 1) or (not d and k == -1):
                                        self.ssind.append('e')
                                    else:
                                        self.ssind.append('b')
                                    d = not d
                                else:
                                    self.sstrans1.append((vh, i))
                                    self.sstrans2.append(
                                        (v1, i1))  # check v or smth

    def write_outs(self, args):
        with open(args.cnct, 'w') as f:
            f.write(''.join(self.cnct))

        with open(args.output, 'w') as f:
            f.write(''.join(self.pdb))

        with open(args.top, 'w') as f:
            f.write(''.join(self.add))

if __name__ == '__main__':

    args_ = get_args()
    cosm = COSM()
    cosm.process(args_)
    cosm.load_sequences(args_.seq, args_.oligs)
    cosm.build_scheme(args_)
    cosm.build_newscheme()
    cosm.assign_ff_particles()
    cosm.process_next_()
    cosm.write_outs(args_)
