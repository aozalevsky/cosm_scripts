#!/usr/bin/env python
# encoding: utf-8

from __future__ import division
#from math import atan2, atan, hypot, pi, sin, cos, acos, sqrt
import numpy as np
import argparse
from random import random, choice
import os
import os.path as osp

#############################################################
#   DOES NOT WORK FOR SINGLE-STRANDED OLIGS (BUT REQUIRES)
#############################################################

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True,
                    action='store', dest='input',
                    help='Coarse-coarse-grained model pdb file')
parser.add_argument('-o', '--output', required=True,
                    action='store', dest='output',
                    help='Atomic model pdb file')
parser.add_argument('-t', '--top', required=True,
                    action='store', dest='top',
                    help='File with ends and crossovers coords')
parser.add_argument('-l', '--lattice', required=True,
                    action='store', dest='lattice',
                    help='Lattice type (honeycomb / square)')
parser.add_argument('-p', '--cgprob', default=None,
                    action='store', dest='cgprob', type=float,
                    help='CG probability')
parser.add_argument('-s', '--seq', default=None,
                    action='store', dest='seqoligs', type=str,
                    help='Oligs csv')

# TER scaffold --> staples

args = parser.parse_args()

if args.lattice in ['h', 'hex', 'hexagonal']:
    SQ = False
elif args.lattice in ['s', 'sq', 'square']:
    SQ = True
else:
    raise Exception("Choose lattice (hexagonal / square)")

if not args.cgprob:
    if args.seqoligs:
        SEQ = True
    else:
        raise Exception("Choose CG ratio or oligs.csv")
else:
    if args.seqoligs:
        raise Exception("Remove CG ratio or oligs.csv")
    else:
        SEQ = False

count = {'T1': 1, 'T2': 2, 'T3': 3, 'T4': 4, 'T5': 5,
         'T6': 6, 'H': 7, 'T7': 7, 'T': 7, 'PT': 7, 'S': 1,
         'TT': 7, 'T1T': 1, 'T2T': 2, 'T3T': 3, 'T4T': 4,
         'T5T': 5, 'T6T': 6, 'T7T': 7, 'O': 1, 'OT': 1, 'N': 1,
         'B1': 1, 'B2': 2, 'B3': 3, 'B4': 4, 'B5': 5, 'B6': 6,
         'B7': 7, 'B1T': 1, 'B2T': 2, 'B3T': 3, 'B4T': 4,
         'B5T': 5, 'B6T': 6, 'B7T': 7, 'B': 7, 'ST': 1}

term = ['T', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'TT', 'T1T', 'T2T',
        'T3T', 'T4T', 'T5T', 'T6T', 'T7T', 'T7',
        'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7',
        'B1T', 'B2T', 'B3T', 'B4T', 'B5T', 'B6T', 'B7T', 'BT', 'B']
if SQ:
    count['PT'] = 8
    count['T'] = 8
    count['B'] = 8
    count['H'] = 8
    count['TT'] = 8

if args.cgprob > 1:
    cgprob = args.cgprob / 100.0
else:
    cgprob = args.cgprob

#rgs()

# -------- load templates -------

print '\tLoad templates...'

# template_path_key
tpk = 'TEMPLATE_PATH'
# template_path
tp = ''

if tpk in os.environ:
    tp = os.environ[tpk]

ncl_a = {}
ncl_c = {}

for template in ['at', 'ta', 'gc', 'cg']:
    res = ''
    with open(osp.join(tp, template + '.pdb'), 'r') as tmp:
        for line in tmp:
            line = line.split()
            if line[0] == 'ATOM':
                ncl_name = line[3] + line[4]
                ncl_a[ncl_name] = ncl_a.get(ncl_name, [])
                ncl_c[ncl_name] = ncl_c.get(ncl_name, [])
                if line[5] != res:
                    i = 0
                    res = line[5]
                if line[2] not in ncl_a[ncl_name]:
                    ncl_a[ncl_name].append(line[2])
                    ncl_c[ncl_name].append([])
                coords = [float(c) for c in line[6:]]
                ncl_c[ncl_name][i] += [coords]
                i += 1

for name in ncl_c:
    if name[1] == 'B':
        for k, a in enumerate(ncl_c[name]):
            ncl_c[name][k] = a[::-1]
    for k, a in enumerate(ncl_c[name]):
        for i, c in enumerate(a):
            c[2] += 3.4 * i

#print len(ncl_a['AA']), ncl_a['AA']
#print len(ncl_c['TB']), ncl_c['TB']

# -------- load pdb ---------

print '\tLoad structure...'


'''
def angle(sbase, nbase):
    """Get angle of base pair for each coord"""
    angle = ((3.508 * sbase) % pi)
#    if sq:
#        angle = ((270 * sbase / 8) % 360) * pi / 180
#    else:
#        angle = ((720 * sbase / 21) % 360) * pi / 180
    print angle
    return angle


def sphere(coords):
    [x, y, z] = coords
    # phi - 2pi, tetha - pi. normally pi/2
    r = sqrt(x ** 2 + y ** 2 + z ** 2)
    phi = atan2(y, x)
    tetha = acos(z / r)
#    print x, y, z, r, phi, tetha
    return [r, phi, tetha]


def cylindr(coords):
    [x, y, z] = coords
    r = hypot(x, y)
    phi = atan2(x, y)
    return [z, r, phi]


def cartc(atom):
    [h, r, phi] = atom
    x = r * sin(phi)
    y = r * cos(phi)
    return [x, y, h]


def carts(atom):
    [r, phi, tetha] = atom
    x = r * sin(tetha) * cos(phi)
    y = r * sin(tetha) * sin(phi)
    z = r * cos(tetha)
    return [x, y, z]
'''


def normal(i, source):
    '''#Angles between normal (i-1, i+1) and common normal y,x=0'''
    if source == scaffold:
        prev = source[i - 1]
        next = source[i + 1]
#        print i-1, i+1, 'xxx'
    else:
        prev = source[i - 1][0]
        next = source[i + 1][0]
#    print '------'
#    print prev
#    print source[i]
#    print next
    norm = np.array([a - b for a, b in zip(next, prev)])
#    print i
#    print prev
#   print next
#    print norm
    norm = norm / np.linalg.norm(norm)
    [nx, ny, nz] = norm
#    print norm
    if ny == 0:
        ny = 0.0001
    x = [1, (- nx - nz) / ny, 1]
#    if x[1] < 0:
#        x[1] = -1 * x[1]
    x = x / np.linalg.norm(x)
    y = np.cross(x, norm)
    y = y / np.linalg.norm(y)
#    A = np.array([norm, x, y]).T
    A = np.array([norm, x, y])
    B = np.linalg.inv(A).T
    [bx, by, bz] = B
    bx = bx / np.linalg.norm(bx)
    by = by / np.linalg.norm(by)
    bz = bz / np.linalg.norm(bz)
    C = np.array([bx, by, bz]).T
#    print C, 'CCCCCC'
#    print C[0][0]**2 + C[1][0]**2 + C[2][0]**2
#    print '--------------'
    O = np.array([[0, 0, -1], [1, 0, 0], [0, 1, 0]])
#    print rnum, A, np.dot(A,O)*-1
    return np.dot(C, O)

# center --> base


def create_base(base, btype, i, norm):
    '''Move and rotate avery atom in the base'''
    res = []
    k = i % 10
#    k = 0
    for n in range(len(ncl_c[btype])):
        natom = ncl_c[btype][n][k]
#        natom = atom
#        natom = cylindr(atom)
#        natom[2] += ang
#        natom = cartc(natom)
#        print atom
        natom = np.array(natom).T
        natom = np.dot(norm, natom)
        natom = natom.tolist()
#        print natom
        natom = [a + b for a, b in zip(natom, base)]
        res.append(natom)
    return res


def print_base(baseatoms, btype, end):
    global anum, rnum
    template = "{0[0]:<6s}{0[1]:>5d}  {0[2]:<4s}{0[3]:<3s} {0[4]:>1s}{0[5]:>4d}    {0[6]:>8.3f}{0[7]:>8.3f}{0[8]:>8.3f}{0[9]:>6.2f}{0[10]:>6.2f}"
    ratoms = ncl_a[btype]
    if btype[1] == 'A':
        chain = 'A'
    else:
        chain = 'B'
    name = 'D' + btype[0]
    if end == 5:
        name += '5'
        ratoms = ratoms[3:]
        baseatoms = baseatoms[3:]
    elif end == 3:
        name += '3'
    for i, atom in enumerate(baseatoms):
        t = tuple(['ATOM'] + [anum] + [ratoms[i]] + [name] + [chain] +
                  [rnum] + atom + [1.00, 0.00])
        pdb.write(template.format(t) + '\n')
        anum += 1
        if anum == 100000:
            anum = 1
    rnum += 1
    if rnum == 10000:
        rnum = 1


def stpath(end):
    [number, add] = end
    s = str(number) + '+' + str(add) + '---'
    tag = 'e3'
    while tag != 'e5':
        f = 0
        if [number, add] in e5:
            s += str(number) + '+' + str(add)
            tag = 'e5'
        else:
            for pair in cross:
                if pair[0] == [number, add]:
                    s += str(number) + '+' + str(add) + '/'
                    [number, add] = pair[1]
                    s += str(number) + '+' + str(add) + '---'
                    f = 1
                    cross.remove(pair)
                elif pair[1] == [number, add]:
                    s += str(number) + '+' + str(add) + '/'
                    [number, add] = pair[0]
                    s += str(number) + '+' + str(add) + '---'
                    f = 1
                    cross.remove(pair)
        if not f:
            add += 1
            if add == atoms[number]:
                add = 0
                if number not in scaf:
                    number += 1
                    if number > num:
                        number = 1
                else:
                    number = scaf[number]
    return stcoords(s)


def fdelta(n1, n2):
    res = []
    [bnum1, add1] = n1
    [bnum2, add2] = n2
    if bnum1 == bnum2:
        for i in range(add1, add2 + 1):
            res.append([scaffold[ascaf[(bnum1, i)]], ascaf[(bnum1, i)]])
    else:
        base = bnum1
        c = atoms[base]
        for i in range(add1, c):
            res.append([scaffold[ascaf[(base, i)]], ascaf[(base, i)]])
        if base in scaf:
            base = scaf[base]
        else:
            base += 1
            if base > num:
                base = 1
        while base != bnum2:
            for i in range(atoms[base]):
                res.append([scaffold[ascaf[(base, i)]], ascaf[(base, i)]])
            if base in scaf:
                base = scaf[base]
            else:
                base += 1
                if base > num:
                    base = 1
        for i in range(add2 + 1):
            res.append([scaffold[ascaf[(base, i)]], ascaf[(base, i)]])
    return res


def stcoords(s):
    res = []
#    print s
    s = s.split('/')
    for part in s:
        part = part.split('---')
        beg = part[0].split('+')
        end = part[1].split('+')
        beg = [int(a) for a in beg]
        end = [int(a) for a in end]
        res += fdelta(beg, end)
    return res


def rd():
    i = random()
    if i > cgprob:
        k = choice(['TA', 'AA'])
    else:
        k = choice(['CA', 'GA'])
    if cur_st:
        k = k[:-1] + 'B'
    return k


bases = {'AA': 'DA', 'AB': 'DA', 'TA': 'DT', 'TB': 'DT', 'CA': 'DC',
         'CB': 'DC', 'GA': 'DG', 'GB': 'DG'}

cart_coords = {}

scaffold = []
ascaf = {}
staples = []

atoms = {}
sccross = []
scafres = None
print '\tCreating base\'s coords...'
with open(args.input, 'r') as f:
    d = 0
    s = 1
    beg = 0
    for i, line in enumerate(f):
        if line[0:4] == 'ATOM':
            name = line[16:21].strip()
            coords = [float(line[30:38]), float(line[38:46]),
                      float(line[46:54])]
            c = count[name]
            num = int(line[6:11])
#            print num, name, beg
            if name in ['S', 'ST', 'NT', 'N']:
                if name[0] == 'S' and beg:
                    beg = 0
                    scaffold.append(begin)
                    atoms[pnum] = 1
                    ascaf[(pnum, 0)] = s - 1
                    s += 1
                scaffold.append(coords)
                atoms[num] = 1
                ascaf[(num, 0)] = s - 1
                s += 1
            elif name in term and not beg:
                begin = coords
                csum = c
                pnum = num
                beg = 1
            elif name in term and beg:
                end = coords
                csum = c
                delta = [(a - b) / csum for a, b in zip(end, begin)]
                for n in range(0, csum):
                    origin = [a + b * n for a, b in zip(begin, delta)]
                    scaffold.append(origin)
                    ascaf[(pnum, n)] = s - 1
                    atoms[pnum] = c
                    s += 1
                scaffold.append(end)
                ascaf[(num, 0)] = s - 1
                s += 1
                atoms[num] = 1
                beg = 0
                scafres = s - 1
            elif name not in ['O', 'OT']:
                end = coords
                delta = [(a - b) / csum for a, b in zip(end, begin)]
                for n in range(csum):
                    origin = [a + b * n for a, b in zip(begin, delta)]
                    scaffold.append(origin)
                    ascaf[(pnum, n)] = s - 1
                    s += 1
                atoms[pnum] = csum
                begin = coords
                csum = c
                pnum = num
            else:
#                if not scafres:
#                    scafres = s - 1
                atoms[num] = 1
                ascaf[(num, 0)] = s - 1
                s += 1
                scaffold.append(coords)
#print len(scaffold)
#print mmm
e3 = []
e5 = []
cross = []
scaf = {}
seq = ''
#print atoms
print '\t"Topology" file reading...'

with open(args.top, 'r') as a:
    for line in a:
        line = line.split()
        if line[0] == 'e5':
            e5.append([int(line[2]), int(line[3])])
        elif line[0] == 'e3':
            e3.append([int(line[2]), int(line[3])])
        elif line[0] == 'cross':
            cross.append([[int(line[2]), int(line[3])], [int(line[4]), int(line[5])]])
        elif line[0] == 'scaf':
            scaf[int(line[2])] = int(line[4])
        elif line[0] == 'seq:':
            tseq = line[1]
if SEQ and not tseq:
    raise Exception("No sequence in topology file")
#            scaf[int(line[4])] = int(line[2])
#print atoms
for end in e3:
    staple = stpath(end)
#    print staple
    staple.reverse()
    staples.append(staple)
# make ends!
#staples.reverse()
#ocnct = {}
#with open(args.connects, 'r') as c:
#    for line in c:
#        line = line.split()
#        ocnct[int(line[0])] = int(line[1])

comp = {'AA': 'TB', 'TA': 'AB', 'CA': 'GB', 'GA': 'CB'}
pdb = open(args.output, 'w')

#angles = {}

NORM = {}
anum = 1
rnum = 1
print '\tCalculate fullatom ...'
cur_st = False
seq = {}
for i, base in enumerate(scaffold[:scafres + 1]):
    end = None
#    print i, 'ooo'
    if i == 0:
        end = 5
        A = normal(1, scaffold)
    elif i == len(scaffold) - 1:
        end = 3
        A = normal(i - 1, scaffold)
    else:
        A = normal(i, scaffold)
    NORM[i] = A
#    print A
#    ang = angle(i, 1)
#    angles[i] = i % 10
#    print i
    if SEQ:
        seq[i] = tseq[i].upper() + 'A'
    else:
        seq[i] = rd()
    nbase = create_base(base, seq[i], i, A)
    print_base(nbase, seq[i], end)

cur_st = True
pdb.write('TER\n')

for staple in staples:
    for i, base in enumerate(staple):
        end = None
        if True:
            if i == 0:
                end = 5
                A = normal(1, staple)
            elif i == len(staple) - 1:
                end = 3
                A = normal(i - 1, staple)
            else:
                A = normal(i, staple)
        if base[1] in NORM:
            A = NORM[base[1]]
#        # why i? check
#        if base[1] in angles:
#            ang = angles[base[1]]
#        else:
#            ang = angles[0]
        if len(seq) - 1 >= base[1]:
            s = comp[seq[base[1]]]
        else:
            s = rd()
        nbase = create_base(base[0], s, base[1], A)
        print_base(nbase, s, end)
    pdb.write('TER\n')
pdb.close()
print '\t\tDone'
