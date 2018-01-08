#!/usr/bin/env python
# encoding: utf-8

from __future__ import division
import numpy as np
import argparse
from random import random, choice
import os


try:
    TEMPLATE_PATH = os.environ['TEMPLATE_PATH']
    print(TEMPLATE_PATH)
except:
    TEMPLATE_PATH = ''


#
#   DOES NOT WORK FOR SINGLE-STRANDED OLIGS (BUT REQUIRES)
#

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


def normal(i, source):
    '''#Angles between normal (i-1, i+1) and common normal y,x=0'''
    if source == scaffold:
        prev = source[i - 1]
        next = source[i + 1]
    else:
        prev = source[i - 1][0]
        next = source[i + 1][0]
    norm = np.array([a - b for a, b in zip(next, prev)])
    if not norm.any():
        return normal(i + 1, source)
    norm = norm / np.linalg.norm(norm)
    [nx, ny, nz] = norm
    if ny == 0:
        # if nx * nz > 0:
        # ny = 0.0001
        # else:
        ny = 0.0001
        # if ny < 0:
        # [nx, ny] = [-nx, -ny]
    x = [1, (- nx - nz) / ny, 1]
    if ny < 0:
        x[1] = x[1] * -1
        x[0] = -1
        x[2] = -1
    x = x / np.linalg.norm(x)
    y = np.cross(x, norm)
    y = y / np.linalg.norm(y)
    A = np.array([norm, x, y])
    B = np.linalg.inv(A).T
    [bx, by, bz] = B
    bx = bx / np.linalg.norm(bx)
    by = by / np.linalg.norm(by)
    bz = bz / np.linalg.norm(bz)
    C = np.array([bx, by, bz]).T
    O = np.array([[0, 0, -1], [1, 0, 0], [0, 1, 0]])
    return np.dot(C, O)

# center --> base


def create_base(base, btype, i, norm):
    '''Move and rotate avery atom in the base'''
    res = []
    k = i % 10
    for n in range(len(ncl_c[btype])):
        natom = ncl_c[btype][n][k]
        natom = np.array(natom).T
        natom = np.dot(norm, natom)
        natom = natom.tolist()
        natom = [a + b for a, b in zip(natom, base)]
        res.append(natom)
    return res


def print_base(baseatoms, btype, end):
    global anum, rnum
    template = (
        "{0[0]:<6s}{0[1]:>5d}  {0[2]:<4s}{0[3]:<3s} {0[4]:>1s}{0[5]:>4d}    "
        "{0[6]:>8.3f}{0[7]:>8.3f}{0[8]:>8.3f}{0[9]:>6.2f}{0[10]:>6.2f}"
        )
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
            if add == atoms[number]:
                add = 0
                if number not in scaf:
                    number += 1
                    if number > len(atoms):
                        number = 1
                else:
                    number = scaf[number]
            else:
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
            res.append([scaffold[ascaf[(bnum1, i)] - 1], ascaf[(bnum1, i)]])
    else:
        base = bnum1
        c = atoms[base]
        for i in range(add1, c):
            res.append([scaffold[ascaf[(base, i)] - 1], ascaf[(base, i)]])
        if base in scaf:
            base = scaf[base]
        else:
            base += 1
            if base > num:
                base = 1
        while base != bnum2:
            for i in range(atoms[base]):
                res.append([scaffold[ascaf[(base, i)] - 1], ascaf[(base, i)]])
            if base in scaf:
                base = scaf[base]
            else:
                base += 1
                if base > num:
                    base = 1
        for i in range(add2 + 1):
            res.append([scaffold[ascaf[(base, i)] - 1], ascaf[(base, i)]])
    return res


def stcoords(s):
    res = []
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


def parse_int(names, coords, nums):
    global rnum_tmp
    if len(names) == 1:
        scaffold.append(coords[0])
        atoms[nums[0]] = 1  # number -> count
        ascaf[(nums[0], 0)] = len(scaffold)  # number+add -> full_number
    else:
        scaffold.append(coords[0])
        ascaf[(nums[0], 0)] = len(scaffold)
        if len(names) == 2:
            if names[0] in ['T', 'TT']:
                c = count[names[1]]
            else:
                c = count[names[0]]
            atoms[nums[0]] = c
            atoms[nums[1]] = 1
            delta = [(a - b) / c for a, b in zip(coords[1], coords[0])]
            for n in range(1, c):
                origin = [a + b * n for a, b in zip(coords[0], delta)]
                scaffold.append(origin)
                ascaf[(nums[0], n)] = len(scaffold)
            scaffold.append(coords[1])
            ascaf[(nums[1], 0)] = len(scaffold)
        else:
            if names[1] in ['H', 'PT']:
                c = count[names[0]]
                atoms[nums[0]] = c
                delta = [(a - b) / c for a, b in zip(coords[1], coords[0])]
                for n in range(1, c):
                    origin = [a + b * n for a, b in zip(coords[0], delta)]
                    scaffold.append(origin)
                    ascaf[(nums[0], n)] = len(scaffold)
                scaffold.append(coords[1])
                ascaf[(nums[1], 0)] = len(scaffold)
            else:
                atoms[nums[0]] = 1
                scaffold.append(coords[1])
                ascaf[(nums[1], 0)] = len(scaffold)
            if len(names) == 3:
                if names[1] in ['H', 'PT']:
                    c = count[names[2]]
                    atoms[nums[1]] = c
                    atoms[nums[2]] = 1
                    delta = [(a - b) / c for a, b in zip(coords[2], coords[1])]
                    for n in range(1, c):
                        origin = [a + b * n for a, b in zip(coords[1], delta)]
                        scaffold.append(origin)
                        ascaf[(nums[1], n)] = len(scaffold)
                    scaffold.append(coords[2])
                    ascaf[(nums[2], 0)] = len(scaffold)
                else:
                    atoms[nums[1]] = 1
                    atoms[nums[2]] = 1
                    scaffold.append(coords[2])
                    ascaf[(nums[2], 0)] = len(scaffold)

            else:
                # print names, nums
                for i in range(2, len(names) - 1):
                    if names[i] not in 'N':
                        c = count[names[i - 1]]
                        atoms[nums[i - 1]] = c
                        delta = [
                            (a - b) / c for a, b in zip(
                                coords[i], coords[i - 1])
                            ]
                        for n in range(1, c):
                            origin = [
                                a + b * n for a, b in zip(coords[i - 1], delta)]
                            scaffold.append(origin)
                            ascaf[(nums[i - 1], n)] = len(scaffold)
                        scaffold.append(coords[i])
                        ascaf[(nums[i], 0)] = len(scaffold)
                    else:
                        atoms[nums[i - 1]] = 1
                        scaffold.append(coords[i])
                        ascaf[(nums[i], 0)] = len(scaffold)
                if names[-2] in ['H', 'PT']:
                    c = count[names[-1]]
                    atoms[nums[-2]] = c
                    atoms[nums[-1]] = 1
                    delta = [
                        (a - b) / c for a, b in zip(coords[-1], coords[-2])]
                    for n in range(1, c):
                        origin = [a + b * n for a, b in zip(coords[-2], delta)]
                        scaffold.append(origin)
                        ascaf[(nums[-2], n)] = len(scaffold)
                    scaffold.append(coords[-1])
                    ascaf[(nums[-1], 0)] = len(scaffold)
                else:
                    atoms[nums[-2]] = 1
                    atoms[nums[-1]] = 1
                    scaffold.append(coords[-1])
                    ascaf[(nums[-1], 0)] = len(scaffold)


# ------------- parameters -----------------

if args.lattice in ['h', 'hex', 'hexagonal']:
    SQ = False
elif args.lattice in ['s', 'sq', 'square']:
    SQ = True
else:
    raise Exception("ADMIN: Wrong lattice type")

if not args.cgprob:
    if args.seqoligs:
        SEQ = True
    else:
        raise Exception("ADMIN: Choose CG ratio or oligs.csv")
else:
    if args.seqoligs:
        raise Exception("ADMIN: Remove CG ratio or oligs.csv")
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

# if args.cgprob > 1:
#    cgprob = args.cgprob / 100.0
# else:
cgprob = args.cgprob

# -------- load templates -------


ncl_a = {}
ncl_c = {}

try:
    for template in ['at', 'ta', 'gc', 'cg']:
        res = ''
        with open(os.path.join(TEMPLATE_PATH, template + '.pdb'), 'r') as tmp:
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

except:
    raise Exception('ADMIN: Error in bases templates files')

# -------- load pdb ---------

# print '\tLoad structure...'

bases = {'AA': 'DA', 'AB': 'DA', 'TA': 'DT', 'TB': 'DT', 'CA': 'DC',
         'CB': 'DC', 'GA': 'DG', 'GB': 'DG'}

cart_coords = {}

scaffold = []
ascaf = {}
staples = []

atoms = {}
sccross = []
scafres = None
# print '\tCreating base\'s coords...'
# try:
if True:
    with open(args.input, 'r') as f:
        d = 0
        s = 1
        beg = 0
        int_name = []
        int_coords = []
        int_num = []
        for i, line in enumerate(f):
            if line[0:4] == 'ATOM':
                name = line[16:21].strip()
                coords = [float(line[30:38]), float(line[38:46]),
                          float(line[46:54])]
                num = int(line[6:11])
                if name in ['S', 'ST']:  # , 'NT', 'N']:
                    parse_int([name], [coords], [num])
                    beg = 0  # really wrong
                else:
                    int_name.append(name)
                    int_coords.append(coords)
                    int_num.append(num)
                    if name in term and not beg:
                        beg = 1
                    elif name in term and beg:
                        parse_int(int_name, int_coords, int_num)
                        beg = 0
                        int_name = []
                        int_coords = []
                        int_num = []

# print ascaf[(9,0)]
# print scaffold[51]
e3 = []
e5 = []
cross = []
scaf = {}
seq = ''
try:
    with open(args.top, 'r') as a:
        for line in a:
            line = line.split()
            if line[0] == 'e5':
                e5.append([int(line[2]), int(line[3])])
            elif line[0] == 'e3':
                e3.append([int(line[2]), int(line[3])])
            elif line[0] == 'cross':
                cross.append(
                    [
                        [int(line[2]), int(line[3])],
                        [int(line[4]), int(line[5])]
                    ])
            elif line[0] == 'scaf':
                scaf[int(line[2])] = int(line[4])
            elif line[0] == 'seq:':
                tseq = line[1]
except:
    raise Exception('ADMIN: PY2. Error in topology file')

if SEQ and not tseq:
    raise Exception(
        "ADMIN: No sequence in topology file, choose -p instead of -s")
# print e3
for end in e3:
    staple = stpath(end)
    staple.reverse()
    staples.append(staple)
# make ends!

comp = {'AA': 'TB', 'TA': 'AB', 'CA': 'GB', 'GA': 'CB'}
pdb = open(args.output, 'w')

NORM = {}
anum = 1
rnum = 1
cur_st = False
seq = {}
for i, base in enumerate(scaffold):
    end = None
    if i == 0:
        end = 5
        A = normal(1, scaffold)
    elif i == len(scaffold) - 1:
        end = 3
        A = normal(i - 1, scaffold)
    else:
        A = normal(i, scaffold)
    NORM[i] = A
    if SEQ:
        seq[i] = tseq[i].upper() + 'A'
    else:
        seq[i] = rd()
    nbase = create_base(base, seq[i], i, A)
    print_base(nbase, seq[i], end)

cur_st = True
pdb.write('TER\n')
for staple in staples:
    # print staple
    for i, base in enumerate(staple):
        end = None
#        if False:
#            if i == 0:
#                end = 5
#                A = normal(1, staple)
#            elif i == len(staple) - 1:
#                end = 3
#                A = normal(i - 1, staple)
#            else:
#                A = normal(i, staple)
        if base[1] - 1 in NORM:
            A = NORM[base[1] - 1]
        else:
            raise Exception('ADMIN: PY2. Single-stranded olig\'s normal error')
        if len(seq) >= base[1]:
            s = comp[seq[base[1] - 1]]
        else:
            raise Exception(
                'ADMIN: PY2. Single-stranded olig\'s sequence error')
            s = rd()
        nbase = create_base(base[0], s, base[1], A)
        print_base(nbase, s, end)
    pdb.write('TER\n')
pdb.close()
