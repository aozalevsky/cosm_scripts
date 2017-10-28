#!/usr/bin/env python
# encoding: utf-8

from __future__ import division
import numpy as np
import argparse
import prody
import re

#
#   DOES NOT WORK FOR SINGLE-STRANDED OLIGS (BUT REQUIRES)
#

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


# In[108]:

CR = [2.004, 2.145, 2.286, 1.0]
SR = [0.33, 0.34, 0.37, 1.0]

LR = [3.68, 3.75, 3.77, 0.1]
CCR = [3.111, 3.204, 3.300, 0.2]

if SQ:
    LR = [4.80, 5.2, 5.23, 0.1]
    CCR = [3.378, 3.464, 3.553, 0.2]


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


def divide_ds(DS, coords):
    DS_ = list()

    for ds in DS:
        ds_ = list()
        x = coords[ds[0]][0]
        y = coords[ds[0]][1]

        for p in ds:
            x_ = coords[p][0]
            y_ = coords[p][1]
            if (x_ == x) and (y_ == y):
                pass
            else:
               DS_.append(ds_)
               ds_ = list()
               x = x_
               y = y_
            ds_.append(p)

    return DS_
                 

        

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
            i, j = int(ts[1]) - 1, int(ts[2]) - 1
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
    nodemap = None
    rmap = None

    cross = None
    junctions = None
    scaf = None
    seq = None

    fname = None
    e3 = None
    e5 = None

    torsions = None

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
                    scaf[int(line[2]) - 1] = int(line[4]) - 1
                elif line[0] == 'seq:':
                    seq = line[1]

    except Exception as e:
        raise Exception('ADMIN: PY2. Error in topology file.' + str(e))

    top.e3 = e3
    top.e5 = e5
    top.cross = cross
    top.scaf = scaf
    top.seq = seq

    # print e5, e3, cross, junctions, cross
    return top


# In[78]:

cgtop = read_cgtop(tfile)


# In[118]:

# Check circle backbone
circle = (0, len(struct) - 1)
rcircle = (len(struct) - 1, 0)
if circle in scross or rcircle in scross:
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

if len(ds) > 0:
    DS.append(ds)

DS = divide_ds(DS, coords)

# In[120]:
# http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
junctions = np.array(sorted([item for sublist in scross for item in sublist]))
# print junctions
LAT = ['H', 'PT', 'T']
lattice = list()
for n in range(len(junctions) - 2):
    i = junctions[n]
    for k in range(n + 1, n + 3):
        j = junctions[k]
        paired = False
        if abs(j - i) <= 2:
            for ds in DS:
                if (i in ds) and (j in ds):
                    if (coords[i][:2] == coords[j][:2]).all():
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

            D2, D3 = None, None

            for ds in DS:
                if ii in ds:
                    D2 = ds
                if jj in ds:
                    D3 = ds

            if D1 and D2 and D3:
                pass
            else:
                continue

            if D1 != D2 != D3:
                pass
            else:
                continue

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
                    if i < 0 or  iin < 0 or jjn < 0:
                        raise IndexError

                    lplane = ((rlabels[i + p] in LAT) and
                              (rlabels[iin] in LAT) and
                              (rlabels[jjn] in LAT)
                              )

                    zline = (
                        (coords[i + p][2]), coords[iin][2], coords[jjn][2])
                    xline = (
                        (coords[i + p][0]), coords[iin][0], coords[jjn][0])
                    yline = (
                        (coords[i + p][1]), coords[iin][1], coords[jjn][1])
                    # print lplane, cplane

                    restr = False
                    if (zline[0] == zline[1] == zline[2]) and lplane:
                        if SQ:
                            if ((xline[0] == xline[1] == xline[2]) or
                                    (yline[0] == yline[1] == yline[2])):
                                restr = True
                        else:
                            restr = True

                    if restr:
                        l = (iin, jjn)
                        rl = (jjn, iin)
                        if l not in lattice and rl not in lattice:
                            lattice.append(l)
                except IndexError as e:
                    pass
                p += 1

cgtop.lattice = sorted(set(lattice))

# Torsions

torsions = list()


for ds_i_ in range(len(DS)):
    ds_i = DS[ds_i_]
    
    if len(ds_i) == 1:
        continue

    for ds_j_ in range(ds_i_):
        ds_j = DS[ds_j_]
        if len(ds_j) == 1:
            continue

        t_s, t_e = None, None

        for s in scross:

            if (s[0] in ds_i and s[1] in ds_j):
            	i_, j_ = s

            elif (s[1] in ds_i and s[0] in ds_j):
            	j_, i_ = s
		
            else:
                continue

            if rlabels[i_] != 'H' or rlabels[j_] != 'H':
                continue

            if not t_s or i_ < t_s[0]:
                t_s = s
            elif not t_e or j_ > t_e[0]:
                t_e = s


	if t_s and t_e:
	        torsions.append((
        	    t_s[1], t_s[0],
	            t_e[0], t_e[1],
        	))

cgtop.torsions = torsions

# In[129]:


def prepare_atom(a, c):
    PDBF = "%-6s%5s %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s  \n"
    coords = a.getCoords() + np.random.normal(0, 0.5, (3, )) / 10
    s = PDBF % (
        'ATOM', c, a.getName(), '',  # altloc empty
        a.getResname(), a.getChid(), a.getResnum(),
        coords[0], coords[1], coords[2],
        1.0, 0.0, '',
        a.getElement()
    )
    return s


def write_fullatom_pdb(fname, top):
    RF = '%d %d 1 %d 1 '

    with open(fname, 'w') as f:
        f.write('\n[ dihedrals ]\n')
        f.write('; torsion angles\n')
        
        for t in top.torsions:
            t = map(lambda x: int(x) + 1, t)
            f.write(' '.join(map(str, t)) + " 1 " + '\n')
        

        N = 1
        f.write('\n[ distance_restraints ]\n')
        f.write('; staple crossovers\n')

        for s in top.scross:
            i, j = s
            f.write(RF % (i + 1, j + 1, N) + ' '.join(map(str, CR)) + '\n')
            N += 1

        f.write('\n; criss cross staple crossovers\n')

        for s in top.ccross:
            i, j = s
            f.write(RF % (i + 1, j + 1, N) + ' '.join(map(str, CCR)) + '\n')
            N += 1

        f.write('\n; scaffold crossovers\n')

        for s in top.scaf.items():
            i, j = s
            f.write(RF % (i + 1, j + 1, N) + ' '.join(map(str, SR)) + '\n')
            N += 1

        f.write('\n; lattice crossovers\n')

        for s in top.lattice:
            i, j = s
            f.write(RF % (i + 1, j + 1, N) + ' '.join(map(str, LR)) + '\n')
            N += 1
        

# In[130]:

fname = args.output
write_fullatom_pdb(fname, cgtop)
