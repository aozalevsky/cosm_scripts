#!/usr/bin/env python
# encoding: utf-8

from __future__ import division
import numpy as np
import argparse
from random import random, choice
import os

#############################################################
#   DOES NOT WORK FOR SINGLE-STRANDED OLIGS (BUT REQUIRES)
#############################################################

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True,
                    action='store', dest='input',
                    help='Coarse-coarse-grained model pdb file')
parser.add_argument('-o', '--output', required=True,
                    action='store', dest='output',
                    help='Gromacs restraints file')
parser.add_argument('-t', '--top', required=True,
                    action='store', dest='top',
                    help='File with ends and crossovers coords')
parser.add_argument('-l', '--lattice', required=True,
                    action='store', dest='lattice',
                    help='Lattice type (honeycomb / square)')

args = parser.parse_args()


# ------------- parameters -----------------

if args.lattice in ['h', 'hex', 'hexagonal']:
    SQ = False
elif args.lattice in ['s', 'sq', 'square']:
    SQ = True
else:
    raise Exception("ADMIN: Wrong lattice type")

# In[626]:

import numpy as np
import prody
import subprocess
import os
import networkx as nx
from transforms3d.axangles import axangle2mat
import copy
import re
from collections import Counter

# In[108]:

SQ = True
LR = [3.68, 3.75, 3.77, 0.1]
SR = [0.33, 0.34, 0.37, 1.0]
CR = [2.004, 2.145, 2.286, 1.0]
CCR = [3.111, 3.204, 3.300, 0.1]

if SQ:
    CCR = [3.378, 3.464, 3.553, 0.1]
    LR = [4.80, 5.2, 5.23, 0.1]



# In[72]:

st = args.input
tfile = args.top

# In[73]:

def is_bck(e):
    i, j = e
    if abs(i - j) == 1:
        return True
    else:
        return False

def remove_period(a, mult=2):
    per = mult * np.pi
    return(a - per * np.floor_divide(a, per))


# In[74]:

def get_cg_struct(path):
    struct = prody.parsePDB(path)
    with open(path, 'r') as f:
        tf = f.readlines()

    totA = len(struct)
    bck, cross = list(), list()
    for s in tf:
        if re.match('CONECT', s):
            ts = s.split()
            # fix pdb to python indexing
            i, j  = int(ts[1]) - 1, int(ts[2]) - 1
            if (-1 < i < totA) and (-1 < j < totA):
                ni = min(i, j)
                nj = max(i, j)
                bond = (ni, nj)
                if is_bck(bond):
                    bck.append(bond)
                else:
                    cross.append(bond)
    bck = sorted(bck)
    cross = sorted(cross)
    bonds = bck + cross
    # print len(struct)
    struct.setBonds(bonds)

    return (struct, bck, cross)

struct, bck, scross = get_cg_struct(st)
coords = struct.getCoords()
# print scross


# In[75]:

rlabels = np.array([r.getResname() for r in struct.iterResidues()])
#  print rlabels


# In[76]:

nodes = np.arange(len(struct))

ssdna = np.where((rlabels == 'S') | (rlabels == 'ST'))[0]
dsdna = np.where((rlabels != 'S') & (rlabels != 'ST'))[0]


# In[77]:

class Topology(object):
    nodemap=None
    rmap = None

    cross = None
    junctions = None
    scaf = None
    seq = None

    fname = None
    e3 = None
    e5 = None


    def __init__(self):
        pass

def read_cgtop(fname):
    # mapping is name of file which ends in _t

    top = Topology()
    top.fname = fname
    e3 = dict()
    e5 = dict()
    scaf = dict()
    cross = list()
    seq = None

    # tfile = '/tmp/hc_linear_t'
    try:
        with open(fname, 'r') as a:
            for line in a:
                line = line.split()
                if line[0] == 'e5':
                    e5[int(line[2]) - 1] = int(line[3])
                elif line[0] == 'e3':
                    e3[int(line[2]) - 1] = int(line[3])
                elif line[0] == 'cross':
                    cross.append(
                        ((int(line[2]) - 1, int(line[3])),
                         (int(line[4]) - 1, int(line[5])))
                    )
                elif line[0] == 'scaf':
                    scaf[int(line[2])] = int(line[4])
                elif line[0] == 'seq:':
                    seq = line[1]

    except Exception, e:
        raise Exception('ADMIN: PY2. Error in topology file.' + str(e))

    top.e3 = e3
    top.e5 = e5
    top.cross = cross
    top.scaf = scaf
    top.seq = seq
    rmap = dict()

    # print e5, e3, cross, junctions, cross
    return top


# In[78]:

cgtop = read_cgtop(tfile)


# In[118]:

# Check circle backbone
circle = (0, len(struct) - 1)
if circle in scross:
    cgtop.scaf.append((0, len(struct)))

ncross = list()
for e in scross:
    i, j = e

    fwc, bkc = False, False
    fwcross = list()
    bkcross = list()

    # Check forward
    if rlabels[j - 1] in ['H', 'PT', 'N']:
        if rlabels[i + 1] in ['H', 'PT', 'N']:
            fwc = True
            fwcross.append((i, j - 1))
            fwcross.append((i + 1, j))
        else:
            fwcross.append((i, j - 1))

    elif rlabels[i + 1] in ['H', 'PT', 'N']:
        fwcross.append((i + 1, j))

    # If found criscross - go to next crossover
    if fwc:
        ncross.extend(fwcross)
        continue

    # Check backward
    if rlabels[j + 1] in ['H', 'PT', 'N']:
        if rlabels[i - 1] in ['H', 'PT', 'N']:
            bkc = True
            bkcross.append((i, j + 1))
            bkcross.append((i - 1, j))
        else:
            bkcross.append((i, j - 1))

    elif rlabels[i - 1] in ['H', 'PT', 'N']:
        bkcross.append((i - 1, j))

    if bkc:
        ncross.extend(bkcross)
        continue

    if len(fwcross) > len(bkcross):
        ncross.extend(fwcross)
    elif len(bkcross) > 0:
        ncross.extend(bkcross)

cgtop.ccross = ncross
cgtop.scross = scross
# print cgtop.cross



# In[119]:

ds_bck = list()
ss_bck = list()

N = 0
DS = list()
ds = list()
while N < len(struct):
    n = N
    if n in ssdna:
        if len(ds) > 0:
            DS.append(ds)

            # print '-' * 20
            # print ds

            ds = list()
    else:
        ds.append(n)
    N += 1


# In[120]:

# http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
junctions = np.array(sorted([item for sublist in scross for item in sublist]))
# print junctions
LAT = ['H', 'PT', 'T']
lattice = list()
for n in range(len(junctions) - 1):
    i = junctions[n]
    j = junctions[n + 1]
    paired = False
    if abs(j - i) <= 2:
        for ds in DS:
            if (i in ds) and (j in ds):
                paired = True
                D1 = ds


        if not paired:
            continue

        for c in scross:
            if (i in c) and (j not in c):
                if i == c[0]:
                    ii = c[1]
                elif i == c[1]:
                    ii = c[0]
            if (i not in c) and (j in c):
                if j == c[0]:
                    jj = c[1]
                elif j == c[1]:
                    jj = c[0]

        for ds in DS:
            if ii in ds:
                D2 = ds
            if jj in ds:
                D3 = ds

        # print i, j, ii, jj
        delta = j - i

        check = True
        R = 3
        p = - R

        while (p < (delta + R + 1)) and check:
            # check backward
            try:

                iin = ii - p
                jjn = jj + delta - p
                # print i, iin, jjn
                lplane = ((rlabels[i + p] in LAT) and
                    (rlabels[iin] in LAT) and
                    (rlabels[jjn] in LAT)
                          )

                cplane = ((coords[i + p][2]), coords[iin][2], coords[jjn][2])
                # print lplane, cplane

                if (cplane[0] == cplane[1] == cplane[2]) and lplane:
                    lattice.append((iin, jjn))
                p += 1
            except IndexError, e:
                pass
cgtop.lattice = lattice
# for e in  sorted(set(lattice)):
#    i, j = e
    # print i + 1, j + 1



# In[129]:

def prepare_atom(a, c):
    PDBF = "%-6s%5s %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s  \n"
    coords = a.getCoords() + np.random.normal(0, 0.5, (3, )) / 10
    s = PDBF % (
        'ATOM', c, a.getName(), '', # altloc empty
        a.getResname(), a.getChid(), a.getResnum(),
        coords[0], coords[1], coords[2],
        1.0, 0.0, '',
        a.getElement()
        )
    return s

def write_fullatom_pdb(fname, top):
    RF = '%d %d 1 %d 1 '

    with open(fname, 'w') as f:
        N = 1
        f.write('[ distance_restraints ]\n')
        f.write('; staple crossovers\n')

        for s in top.scross:
            i, j = s
            f.write(RF % (i + 1, j + 1, N) + ' '.join(map(str, CR)) + '\n')
            N += 1

        f.write('; criss cross staple crossovers\n')

        for s in top.ccross:
            i, j = s
            f.write(RF % (i + 1, j + 1, N) + ' '.join(map(str, CCR)) + '\n')
            N += 1

        f.write('; scaffold crossovers\n')

        for s in top.scaf.items():
            i, j = s
            f.write(RF % (i + 1, j + 1, N) + ' '.join(map(str, SR)) + '\n')
            N += 1

        f.write('; lattice crossovers\n')

        for s in top.lattice:
            i, j = s
            f.write(RF % (i + 1, j + 1, N) + ' '.join(map(str, LR)) + '\n')
            N += 1



# In[130]:

fname = args.output
write_fullatom_pdb(fname, cgtop)


# In[ ]:



