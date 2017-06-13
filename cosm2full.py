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
        'B1' : 1, 'B2': 2, 'B3': 3, 'B4': 4, 'B5': 5, 'B6': 6,
        'B7' : 7, 'B1T': 1, 'B2T' : 2, 'B3T' : 3, 'B4T' : 4,
        'B5T' : 5, 'B6T': 6, 'B7T' : 7, 'B': 7, 'ST': 1}
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

#if args.cgprob > 1:
#    cgprob = args.cgprob / 100.0
#else:
cgprob = args.cgprob

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

import sys
sys.path.append('/home/www-data/web2py/applications/cosm/private/origami_scripts')
import linpol


# In[629]:

RADP = np.deg2rad(360.0 / 10.5)
count = {'T1': 1, 'T2': 2, 'T3': 3, 'T4': 4, 'T5': 5,
         'T6': 6, 'H': 7, 'T7': 7, 'T': 7, 'PT': 7, 'S': 1,
        'TT': 7, 'T1T': 1, 'T2T': 2, 'T3T': 3, 'T4T': 4,
        'T5T': 5, 'T6T': 6, 'T7T': 7, 'O': 1, 'OT': 1, 'N': 1,
        'B1' : 1, 'B2': 2, 'B3': 3, 'B4': 4, 'B5': 5, 'B6': 6,
        'B7' : 7, 'B1T': 1, 'B2T' : 2, 'B3T' : 3, 'B4T' : 4,
        'B5T' : 5, 'B6T': 6, 'B7T' : 7, 'B': 7, 'ST': 1}
if SQ:
    RADP = np.deg2rad(360.0 / 10.67)
    count['PT'] = 8
    count['T'] = 8
    count['B'] = 8
    count['H'] = 8
    count['TT'] = 8


# In[630]:

st = args.input

tfile = args.top

# In[631]:

def is_bck(e):
    i, j = e
    if abs(i - j) == 1:
        return True
    else:
        return False

def remove_period(a, mult=2):
    per = mult * np.pi
    return(a - per * np.floor_divide(a, per))


# In[632]:

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
    print len(struct)
    struct.setBonds(bonds)

    return (struct, bck, cross)

struct, bck, cross = get_cg_struct(st)
cross
# http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
junctions = np.array(sorted([item for sublist in cross for item in sublist]))
print junctions
print cross
fcoords = struct.getCoords()
fcoords += np.random.normal(0, 0.5, fcoords.shape) / 1000


# In[633]:

rlabels = np.array([r.getResname() for r in struct.iterResidues()])
print rlabels


# In[634]:

rcount = np.array([x for x in map(lambda x: count[x], rlabels)])
print sum(rcount)


# In[635]:

rcount = np.zeros(rlabels.shape, dtype=np.int)
rcount[0] = count[rlabels[0]]
for i in range(1, len(rcount)):
    if rlabels[i] in ['T', 'PT'] and rlabels[i - 1] == 'H':
        print i, rlabels[i]
        rcount[i] = 0
    elif rlabels[i][0] == 'T' and rlabels[i - 1] == 'PT':
        print i, rlabels[i]
        rcount[i] = 1
        rcount[i - 1] = count[rlabels[i]]
    else:
        rcount[i] = count[rlabels[i]]
NR = sum(rcount)
print NR


# In[636]:

get_ipython().magic('matplotlib nbagg')
X = np.array(map(lambda x: x[0], fcoords))
Y = np.array(map(lambda x: x[1], fcoords))
Z = np.array(map(lambda x: x[2], fcoords))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(X, Y, Z)


# In[748]:

nodes = np.arange(len(struct))

ssdna = np.where((rlabels == 'S') | (rlabels == 'ST'))[0]
dsdna = np.where((rlabels != 'S') & (rlabels != 'ST'))[0]

topologyd = dict()
ia = 0
ib = sum(rcount[dsdna]) - 1
for i in nodes:
    td = list()
    if i in ssdna:
        for t in range(rcount[i]):
            td.append((ia, ))
            ia += 1

    else:
        for t in range(rcount[i]):
            td.append((ia, ib,))
            ia += 1
            ib -= 1
    topologyd[i] = td


# In[750]:

ds_bck = list()
ss_bck = list()
nfcoords = list()
N = 0
DS = list()
ds = list()
while N < len(topologyd):
    n = N
    l = topologyd[n]
    if len(l[0]) == 1:

        if len(ds) > 0:
            DS.append(ds)
            dss = ds[0]
            dse = ds[-1]
            tcoords = linpol.interparc(
                sum(rcount[ds]),
                    X[ds],
                    Y[ds],
                    Z[ds],
                    )
            # print '-' * 20
            # print ds
            # print tcoords
            nfcoords.extend(tcoords)
            ds = list()
        nfcoords.append(fcoords[n])
    else:
        ds.append(n)
    N += 1
nfcoords = np.array(nfcoords)


# In[640]:

get_ipython().magic('matplotlib nbagg')
nX = np.array(map(lambda x: x[0], nfcoords))
nY = np.array(map(lambda x: x[1], nfcoords))
nZ = np.array(map(lambda x: x[2], nfcoords))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(nX, nY, nZ)


# In[641]:

class Particle(object):
    anum = None
    bnum = None
    angle = None
    coord = None
    single = False
    struct = None
    direction = None

    def __init__(self):
        pass
    def toggle_single(self):
        self.single = not self.single
    def set_ss(self):
        self.single = True
    def set_ds(self):
        self.single = False
    def is_single(self):
        return self.single
    def getCoord(self):
        return self.coord



# In[642]:

fpart = np.ndarray((NR,), dtype=np.object)
for i in topologyd.items():
    ii, l = i
    for j in range(len(l)):
        jj = l[j]
        tp = Particle()
        anum = jj[0]
        tp.anum = anum

        if len(jj) == 2:
            bnum = jj[1]
            tp.bnum = bnum

        tp.coord = nfcoords[anum]
        try:
            d = nfcoords[anum + 1] - nfcoords[anum]
        except IndexError:
            d = nfcoords[anum] - nfcoords[anum - 1]

        tp.direction = d

        if len(jj) == 1:
            tp.set_ss()
        else:
            tp.set_ds()

        fpart[anum] = tp


# In[741]:

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


# In[742]:

def cgtop2ftop(cgtop):
    fe5 = list()
    fe3 = list()
    fcross = dict()


    for i in cgtop.e5.items():
        k, v = i
        fe5.append(cgtop.nodemap[k][v][1])

    for i in cgtop.e3.items():
        k, v = i
        fe3.append(cgtop.nodemap[k][v][1])

    # print fe5, fe3

    for c in cgtop.cross:
        i, j = c

        k, v = i
        ii = cgtop.nodemap[k][v][1]

        k, v = j
        jj = cgtop.nodemap[k][v][1]
        fcross[ii] = jj
        fcross[jj] = ii

    top = Topology()
    top.e5 = fe5
    top.e3 = fe3
    top.cross = copy.copy(fcross)
    top.junctions = set(fcross.keys())

    staples = list()

    for s in fe5:
        staple = list()
        e = s
        staple.append(e)
        e += 1
        while e not in fe3:
            staple.append(e)
            if e in fcross:
                e = fcross.pop(e)
                fcross.pop(e)
            else:
                e += 1
        staple.append(e)

        # print staple
        staples.append(staple)

    top.staples = staples

    rmap = dict()
    rmap['A2B'] = dict()
    rmap['B2A'] = dict()
    rmap['B2N'] = dict()

    for n in cgtop.nodemap.values():
        for b in n:
            if len(b) == 2:
                A, B = b
                rmap['A2B'][A] = B
                rmap['B2A'][B] = A

    top.rmap = rmap

    return top


# In[743]:

cgtop = read_cgtop(tfile)
cgtop.nodemap = topologyd
ftop = cgtop2ftop(cgtop)

if SEQ and not cgtop.seq:
    raise Exception("ADMIN: No sequence in topology file, choose -p instead of -s")

# In[732]:

ds_bck = list()
ss_bck = list()
for p in fpart:
    if p.is_single():
        ss_bck.append((p.anum, p.anum + 1))
    else:
        ds_bck.append((p.anum, p.anum + 1))

dsdna = np.arange(NR, dtype=np.int)[np.nonzero(map(lambda x: not x.is_single(), fpart))]


# In[747]:

tangles = dict()
junc  = set(ftop.junctions)
junc = set(map(lambda x: ftop.rmap['B2A'][x], junc))

cross = map(lambda x: (
        ftop.rmap['B2A'][x[0]],
        ftop.rmap['B2A'][x[1]]),
        ftop.cross.items())

G = nx.Graph()
G.add_nodes_from(dsdna)
G.add_edges_from(ds_bck + cross)

for g in nx.connected_component_subgraphs(G):

    ta = dict()

    gjunct = sorted(set(g.nodes()) & junc)
    # print gjunct

    # for e in nx.bfs_edges(g):
    i = gjunct[0]
    j = ftop.rmap['A2B'][i]
    j = ftop.cross[j]
    j = ftop.rmap['B2A'][j]


    # Find angle of rotation to align horizon with X axis
    v = fpart[i + 1].coord - fpart[i].coord
    v /= np.linalg.norm(v)

    vc1 = fpart[j].coord - fpart[i].coord
    vc1 /= np.linalg.norm(vc1)

    vc2 = fpart[j - 1].coord - fpart[i + 1].coord
    vc2 /= np.linalg.norm(vc2)

    vc = (vc1 + vc2) / 2.0


    vx = np.array([1.0, 0.0, 0.0])

    hza = np.deg2rad(1)
    Mhz = axangle2mat(v, hza)
    delta = 9999999

    for a in range(0, 360, 1):
        tdelta = np.linalg.norm(vx + vc)

        if tdelta < delta:
            ia = a
            delta = tdelta

        vx = np.dot(Mhz, vx)

    ia = np.deg2rad(ia) # - 0.5 * RADP # - RADP * 1.75
    # print i, ia, delta
    ta[i] = ia # + np.pi + RADP


    for e in nx.edge_dfs(g, gjunct[0]):
        i, j = e

        ia = ta[i]
        if is_bck(e):
            if j > i:
                ja = ia + RADP
            else:
                ja = ia - RADP
        else:
            ja = ia - np.pi  # + np.pi

        ta[j] = ja

    tangles.update(ta)
    # print 'WIN'

for e in ss_bck:
    i, j = e
    if i in tangles:
        if j in dsdna:
            continue
        ia = tangles[i]
        ja = ia + RADP
        tangles[j] = ja
    else:
        jt = j
        while jt not in tangles and jt < len(nodes):
            jt += 1
        ja = tangles[jt]

        ia = ja - (jt - i) * RADP
        tangles[i] = ia


angles = np.array(tangles.values())
# angles
# print angles
# angles = {k:remove_period(v) for k, v in angles.items()}


# In[720]:

def gen_duplex(l, s=None, single=False):
    if not s:
        s = np.random.choice(['A', 'T', 'G', 'C'], l, p=[0.25, 0.25, 0.25, 0.25])
    tf = '/tmp/tmpdna.pdb'
    # call = ['/home/domain/silwer/prog/x3dna-v2.1/bin/fiber']
    # call.append('-b') # bform dna
    # call.append('-seq=%s' % ''.join(s)) # sequence
    call = ['/tmp/build_dna.sh']
    call.append('%s' % ''.join(s)) # sequence
    call.append(tf) #file
    # environ = os.environ
    # environ['X3DNA']= '/home/domain/silwer/prog/x3dna-v2.1'
    # subprocess.call(' '.join(call), env=environ, shell=True)
    subprocess.call(call)

    struct = prody.parsePDB(tf)
    if single:
        struct = struct.select('chain A')
        prody.writePDB(tf, struct)
        struct = prody.parsePDB(tf)

    return struct


# In[721]:

tts = gen_duplex(20)


# In[722]:

def move_duplex(o, v, a, struct, single=False):

    vz = np.array([0.0, 0.0, 1.0])
    vx = np.array([1.0, 0.0, 0.0])

    nr = struct.numResidues()
    f = 'resnum %d and chain A' % (1)
    l = 'resnum %d and chain B' % (nr)
    f += ' or ' + l

    # Put A1-BN pair to 0.0 with parallel shift
    ts = struct.select(f)
    ttc = (np.sum(ts.getCoords(), axis=0) / len(ts)) * (-1.0)
    tc = struct.getCoords() + ttc
    struct.setCoords(tc)

    # Calculate axis and angle or rotation to origami vector
    # axis of rotations
    tv = np.cross(vz, v)
    # angle
    avz = np.arccos(np.dot(vz, v) / (np.linalg.norm(vz) * np.linalg.norm(v)))
    # matrix for rotation
    Md = axangle2mat(tv, avz)


    # Fix horizon level for every duplex
    # Use A1 C1' - BN C1' vector as horizon level
    ts = struct.select(f)
    hzv = ts.select("name C1'").getCoords()
    dhz = np.dot(Md, (hzv[1] - hzv[0]))
    dhz /= np.linalg.norm(dhz)

    # Find angle of rotation to align horizon with X axis
    hza = np.deg2rad(1)
    Mhz = axangle2mat(v, hza)
    delta = 999999
    for i in range(0, 360, 1):
        tdelta = np.linalg.norm(vx + dhz)
        if tdelta < delta:
            hza = i
            delta = tdelta
        dhz = np.dot(Mhz, dhz)

    hza = np.deg2rad(hza)
    Mhz = axangle2mat(v, hza)

    # Matrix for final duplex rotation
    Ma = axangle2mat(v, a)

    for i in range(len(tc)):
        c = tc[i]
        c = np.dot(Md, c)
        c = np.dot(Mhz, c)
        c = np.dot(Ma, c)
        tc[i] = c

    tc += o # fix coordinates to origin with parallel shift
    struct.setCoords(tc)

    if single:
        tf = '/tmp/tmpdna.pdb'
        struct = struct.select('chain A')
        prody.writePDB(tf, struct)
        struct = prody.parsePDB(tf)

    return struct


# In[723]:

# tts = gen_duplex(1, single=False)
tts = gen_duplex(10)


# In[724]:

ntts = move_duplex(np.array([0, 0, 0]), np.array([-0.1, -0.1, -1]), 90.0, tts)
prody.writePDB('/tmp/rotdna.pdb', ntts)


# In[725]:

ch = dict()
for i in ['A', 'T', 'G', 'C']:
    ch[i] = gen_duplex(1, s=i)


# In[726]:

for i in range(NR):
    tp = fpart[i]
    tp.angle = angles[i]
    s = np.random.choice(['A', 'T', 'G', 'C'], 1, p=[0.25, 0.25, 0.25, 0.25])[0]
    td = copy.copy(ch[str(s)])
    td = move_duplex(tp.coord, tp.direction, tp.angle, td, single=tp.is_single())
    tp.struct = td


# In[727]:

def merge_particles(parray):
    NR = len(np.nonzero(map(lambda x: not x.is_single(), parray))[0])
    ca = 1
    cb = NR
    print cb

    tff = parray[0].struct.copy()

    for r in tff.iterResidues():
        if r.getChid() == 'A':
            r.setResnum(ca)
            ca += 1
        else:
            r.setResnum(cb)
            cb -= 1

    for i in range(1, len(parray)):

        tt = parray[i].struct.copy()
        for r in tt.iterResidues():
            if r.getChid() == 'A':
                r.setResnum(ca)
                ca += 1
            else:
                r.setResnum(cb)
                cb -= 1

        tff += tt
    return tff

def prepare_scaf(struct):

    tf = '/tmp/tmpdna.pdb'
    struct = struct.select('chain A and (not (resnum 1 and name P OP1 OP2))')
    prody.writePDB(tf, struct)
    struct = prody.parsePDB(tf)

    return struct


# In[728]:

mparts = merge_particles(fpart)
chA = prepare_scaf(mparts)


# In[729]:

def prepare_staps(chB, staples):

    sp = list()
    tc = 0
    ttc = list()
    for ns in staples:
        #print '-' * 20
        #print ns

        tc += len(ns)
        ttc.extend(ns)
        struct = chB.select('(chain B and resnum %s) and not (resnum %d and name P OP1 OP2)' % (
                ' '.join(map(lambda x: str(x + 1), ns)), ns[-1] + 1))

        tf = '/tmp/tmpdna.pdb'
        prody.writePDB(tf, struct)
        stap = prody.parsePDB(tf)

        for r in stap.iterResidues():
            r.setResnum(len(ns) - ns.index(r.getResnum() - 1))

        sp.append(stap)
    #print tc
    #print len(set(ttc))
    #CC = Counter(ttc)
    #print CC
    #for i in range(1, 406):
    #    if i not in CC:
    #        print "XXX", i
    #print sum(map(lambda x: len(x), sp))
    return sp

staps = prepare_staps(mparts, ftop.staples)


# In[744]:

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

def write_fullatom_pdb(scaf, stap, fname):
    with open(fname, 'w') as f:
        c = 1
        for a in scaf.iterAtoms():
            f.write(prepare_atom(a, c))
            c += 1
        f.write('TER\n')

        resi = 0
        for s in stap:
            ts = len(set(s.getResnums()))
            # print ts
            for r in s.iterResidues():
                r.setResnum(resi + r.getResnum())

            atoms = [a for a in s.iterAtoms()]
            atoms = sorted(atoms, key=lambda a: a.getResnum())
            for a in atoms:
                f.write(prepare_atom(a, c))
                c += 1
            resi += ts
            f.write('TER\n')
        # print resi


# In[745]:

write_fullatom_pdb(chA, staps, args.output)
