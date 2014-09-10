#!/usr/bin/env python
# encoding: utf-8

import __main__
__main__.pymol_argv = ['pymol', '-qc']
import pymol
import json
import argparse
import string

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True,
                    action='store', dest='input',
                    help='input cadnano json file')
parser.add_argument('-o', '--output', required=True,
                    action='store', dest='output',
                    help='output pdb file')
parser.add_argument('-r' '--restr', required=True,
                    action='store', dest='cnct',
                    help='output distance/angle restraints')
parser.add_argument('-t' '--top', required=True,
                    action='store', dest='top',
                    help='output file for ccg --> cg / fullatom convertation')
parser.add_argument('-l' '--lattice', required=True,
                    action='store', dest='lattice',
                    help='lattice type (hexagonal / square)')
parser.add_argument('--seq',
                    action='store', dest='seq',
                    help='scaffold sequence')
parser.add_argument('--oligs',
                    action='store', dest='oligs',
                    help='staple sequences exported csv')
args = parser.parse_args()


### do not 3-5 restr if additional oligs
### lattice crossovers - not only near
### other constants for olig-ssnear (hexcl)


def honeycomb(row, col):
    """Cadnano honeycomb lattice into cartesian"""
    # k = distance between helixes
    k = 23.0 / 20.0
    x = col * 17.32 * k
    if (row % 2) != 0:
        y = (row * 3 + 1 - col % 2) * 10 * k
    else:
        y = (row * 3 + col % 2) * 10 * k
    return x, y


def square(row, col):
    """Cadnano square lattice into cartesian"""
    x = col * 23.0
    y = row * 23.0
    return x, y


def revers(vh, scaf):
    """Directoin of scaffold or staple chain (1/-1)"""
    rev = (Rows[vh] + Cols[vh]) % 2
    if rev:
        k = -1
    else:
        k = 1
    if scaf:
        k *= -1
    return k


def get_scaffold_path(end):
#    print end
    #global k
    """Get coords from 5' (end) until find the 3'-end; for scaffold"""
    k = revers(end[0], 0)
#    rev = (Rows[end[0]] + Cols[end[0]]) % 2
    if k == 1:
        path = Path['scaf'][end[0]][end[1] + 1:]
    else:
        path = Path['scaf'][end[0]][:end[1]]
        path.reverse()
    cross_base = next(i for i, str in enumerate(path) if str != end[0])
    cross_vh = path[cross_base]
#    if not rev:
#        k = 1
#    else:
#        k = -1
    cross_base = end[1] + k * cross_base + k
    if cross_vh == 'e3':
        res = ['end', cross_base]
    else:
        for tr in trans_sc:
            if (tr[0] == [end[0], cross_vh]) and (cross_base == tr[1][0]):
                res = [cross_vh, tr[1][1]]
            if (tr[0] == [cross_vh, end[0]]) and (cross_base == tr[1][1]):
                res = [cross_vh, tr[1][0]]
    return [[end[0], cross_base], res]


def nextbase(vh, base, seq):
    """ Get next base after """
    ### works two times in sequnce search:
    ### 1st iteration for next sequence negin point
    ### 2nd iterantion for sequence checking
    global last1, last2
    if seq:
        lastc = last2
    else:
        lastc = last1
#    print vh, base, 'searchnext', lastc
    k = revers(vh, 0)
    if not lastc:
        for pair in trans_sc:
            if pair[0][0] == vh and pair[1][0] == base:
                if seq:
                    last2 = True
                else:
                    last1 = True
                return [pair[0][1], pair[1][1]]
    if seq:
        last2 = False
    else:
        last1 = False
    return [vh, base + k]


def checkseq(vh, base):
    """ Check if sequence and oligs fit on this step """
    global last2
    last2 = last1
    beg = [vh, base]
    k = revers(vh, 0)
    i = 0
    while [vh, base] != beg or not i:
        if allPath[vh][base] in ['A', 'T', 'G', 'C']:
            if seq[i] != compl[allPath[vh][base]]:
                print 'fail'
                return False
        [vh, base] = nextbase(vh, base, True)
        i += 1
    return beg


def findnextnot(end, pat):
    global k
    if k == 1:
        path = allPath[end[0]][end[1]:]
    else:
        path = allPath[end[0]][:end[1]]
        path.reverse()
    if path:
        i = 0
        cross_base = None
        while i < len(path) and not cross_base:
#            print path[i]
            if path[i] != pat:
                cross_base = i
            i += 1
#    cross_vh = path[cross_base]
        if cross_base:
            cross_base = end[1] + k * cross_base
            if k == -1:
                cross_base -= 1
            return [end, [end[0], cross_base]]
        else: return None
    else:
        return None


def pathandtype(end, s):
    """Get next duplex end or vh (chain) change"""
    global k
    if s == 'change':
        new_chain = get_scaffold_path([end[0], end[1] - k])
    else:
        new_chain = get_scaffold_path(end)
    new_change = findnextnot(end, allPath[end[0]][end[1]])
    if new_change:
        if new_chain[0][1] * k <= (new_change[1][1] - k) * k:
            new_end = new_chain[1]
            s = 'chain'
            l = new_chain[0][1]
        else:
            new_end = new_change[1]
            s = 'change'
            l = new_change[1][1] - k
    else:
        new_end = new_chain[1]
        s = 'chain'
        l = new_chain[0][1]
    return [([allPath[end[0]][end[1]], end[0], [end[1], l]], new_end), s]


def prconnect(base1, base2):
    """Print crossover reatraints"""
    global l
    cnct.write(str(base2) + '\t' + str(base1) + '\t1\t' + str(l) + '\t1\t1.7\t1.85\t1.9\t1.0\n')
    l += 1


def findconnects(vh, j, number, length, end):
    """Find crossover restraints"""
    global connects, needed
#    print connst, needst
    if end == 's':
        return None
    [x, y] = [0, 0]
    k = revers(vh, 0)
    seq = [j]
    seq += range(j + k, j + k * length, k)
    ind = []
    for i in seq:
        dif = (i - j) * k
        if Path['stap'][vh][i] not in [vh, 'e5', 'e3', '-']:
            if [vh, i] in needst:
                n = needst.index([vh, i])
                if i == j:
                    needst[n] = (number, 0)
                    add.write('cross ; ' + str(number) + ' 0 ' + str(connst[n][0]) + ' ' + str(connst[n][1]) + '\n')
                elif i % c == c - 1:
                    needst[n] = (number, dif)
                    add.write('cross ; ' + str(number) + ' ' + str(length - 1) + ' '  + str(connst[n][0]) + ' ' + str(connst[n][1]) + '\n')
                else:
                    needst[n] = (number + dif, 0)
                    add.write('cross ; ' + str(number + dif) + ' 0 ' + str(connst[n][0]) + ' ' + str(connst[n][1]) + '\n')
            if [vh, i] in needed:
                n = needed.index([vh, i])
                if i == j:
                    needed[n] = number
                    prconnect(connects[n], number)
                elif i % c == c - 1:
                    if k == 1:
                        needed[n] = number + 1
                        prconnect(connects[n], number + 1)
                    else:
                        needed[n] = number
                        prconnect(connects[n], number)
                else:
                    ind.append(dif)
                    prconnect(connects[n], number + dif) # and here the function?
#                    add.write('cross; ' + str(number + dif) + '0 ' + str(connects[n]) + ' 0\n')
                    needed[n] = number + dif
#                    add.write('cross; ' + str(number) + ' ' + str(c - 1) + ' '  + str(connst[n][0]) + ' ' + str(connst[n][1]) + '\n')
                #add.write
    for i in seq:
        dif = (i - j) * k
        if Path['stap'][vh][i] not in [vh, 'e5', 'e3', '-']:
            if [vh, i] not in needed:
                l = Path['stap'][vh][i]
                for pair in trans:
                    if pair[0] == [vh, l] and pair[1][0] == i:
                        needst.append([pair[0][1], pair[1][1]])
                        if i == j:
                            connst.append((number, 0))
                            if number not in connects + needed:
                                connects.append(number)
                                needed.append([pair[0][1], pair[1][1]])
                        elif i % c == c - 1:
                            connst.append((number, dif))
                            if [vh, i + 1] not in needed:
                                if k == 1:
                                    if number + 1 not in connects + needed:
                                        connects.append(number + 1)
                                        needed.append([pair[0][1], pair[1][1]])
                                else:
                                    if number not in connects + needed:
                                        connects.append(number)
                                        needed.append([pair[0][1], pair[1][1]])
                        else:
                            connst.append((number + dif, 0))
                            if dif not in ind: ind.append(dif)
                            connects.append(number + dif)
                            needed.append([pair[0][1], pair[1][1]])
#                            connst.append((number - 1, c - 1))
#        elif Path['stap'][vh][i] in ['e5', 'e3']:
#            en = Path['stap'][vh][i]
#            if dif != -1:
#                add.write(en + ' ; ' + str(a_sc) + ' ' + str(dif) + '\n')
#            else:
#                add.write(en + ' ; ' + str(a_sc - 1) + ' ' + str(c - 1) + '\n')

    if ind:
        return ind
    return None
        #### if crossover at 6/7 in another chain?
        #### if c-1 - N?

# e5-e3 restr only if circular!


def staplesends(vh, i, number, length):
    """Find staple ends through path"""
    global needst, connst
    rev = (Rows[vh] + Cols[vh]) % 2
    if not rev:
        k = 1
    else:
        k = -1
    l = 0
    while l < length:
        a = Path['stap'][vh][i]
        if a == 'e5':
            add.write('e5 ; ' + str(number) + ' ' + str(l) + '\n')
        elif a == 'e3':
            add.write('e3 ; ' + str(number) + ' ' + str(l) + '\n')
        i += k
        l += 1


def scaffoldcross(vh, i, number):
    """Find and print scaffold crossovers"""
    global scconn, scneed
    rev = (Rows[vh] + Cols[vh]) % 2
    if not rev:
        k = 1
    else:
        k = -1
    if (vh, i) in scneed:
        n = scneed.index((vh, i))
        scneed[n] = number
    elif (vh, i) in scconn:
        n = scconn.index((vh, i))
        scconn[n] = number
    else:
        for pair in trans_sc:
            if pair[0][0] == vh and pair[1][0] == i:
                if Path['scaf'][vh][i + k] not in [vh, '-']:
                    scneed.append((vh, i + k))
                    scconn.append(number)
                elif Path['scaf'][vh][i - k] not in [vh, '-']:
                    scconn.append((vh, i - k))
                    scneed.append(number)


def printatom(aname, tname, nlname, ox, oy, z, chain, length):
    """Print atom in pdb"""
    global a_sc
    template = "{0[0]:<6s}{0[1]:>5d}  {0[2]:<4s}{0[3]:<3s} {0[4]:>1s}{0[5]:>4d}    {0[6]:>8.3f}{0[7]:>8.3f}{0[8]:>8.3f}{0[9]:>6.2f}{0[10]:>6.2f}            {0[11]:<2s}"
    t = tuple(['ATOM'] + [a_sc] + [aname] + [tname] + [nlname] + [a_sc] + [ox, oy]
              + [z * 3.5] + [1.00, 0.00] + [aname])
    pdb.write(template.format(t) + '\n')
    staplesends(chain, z, a_sc, length)
    scaffoldcross(chain, z, a_sc)
    a_sc += 1


def createatom(z, chain, tname, length, end):
    """Print particles in the interval -- calculate coords etc"""
    global a_sc, dif
    [row, col] = [Rows[chain], Cols[chain]]
#    pirint chain, z
    if not SQ:
        [ox, oy] = list(honeycomb(row, col))
    else:
        [ox, oy] = list(square(row, col))
    if a_sc == 1:
        dif = None
    aname = 'D'
    atomsh.append((chain, z))
    anames.append(tname)
    k = revers(chain, 0)
    if (chain, z) in ssnear:
        n = ssnear.index((chain, z))
        ssnear[n] = a_sc
    elif (chain, z - 1) in ssnear:
        n = ssnear.index((chain, z - 1))
        ssnear[n] = a_sc
    if (chain, z) in ssoligs:
        n = ssoligs.index((chain, z))
        if ssnear[n]:
            add.write('scaf ; ' + str(ssnear[n]) + ' ; ' + str(a_sc) + '\n')
    nlname = nls[nl]
    tlname = tname
    if tname == 'N':
        printatom(aname, 'N', nlname, ox, oy, z, chain, 1)  # print untyp
#        dif = findconnects(chain, z, a_sc, length, end) ### problems with already N atoms + connects. but they don't have any connects!
    else:
        if dif and tname not in term:
            tname = 'PT'
        if dif and tname in term and tname[0] == 'D':
            tname = 'T'
        dif = findconnects(chain, z, a_sc, length, end)
        if dif:
            if tname not in term: tname = 'PT'
            if tname in term and tname[0] == 'D': tname = 'T'
            printatom(aname, tname, nlname, ox, oy, z, chain, 1)  # print untyp
            for i in range(1, length):
                printatom(aname, 'N', nlname, ox, oy, z + i * k, chain, 1)  # print untyp
                atomsh.append((chain, z + i * k))
                anames.append('N')
        else:
            printatom(aname, tname, nlname, ox, oy, z, chain, length)  # print this atom. if untyp - PT, if norm - norm


# ------------- .json analysis ----------------

file = open(args.input, 'r')
lines = file.readlines()
stringl = ""
for line in lines:
    stringl += line
obj = json.loads(stringl)
strands = obj["vstrands"]
name = obj["name"]

# create dictionaries (keyed by virtual helix #) of
# row/col, scaf array, stap array

vhToScaf = {}
vhToStap = {}
vhNums = []
Rows = {}
Cols = {}
trans = []
trans_sc = []
end5scaf = None
end3scaf = None
for strand in strands:
    num = strand["num"]
    vhNums.append(num)
    Rows[num] = strand["row"]
    Cols[num] = strand["col"]
    scaf = strand["scaf"]
    stap = strand["stap"]
    vhToScaf[num] = scaf
    vhToStap[num] = stap

scaf5 = []
scaf3 = []
coords = []

Path = {'scaf': {}, 'stap': {}}
ends5 = []
ends3 = []
allPath = {}

for vh in vhNums:
    scafPath = []
    stapPath = []
    scaf = vhToScaf[vh]

    # scaffold path

    for i in range(len(scaf)):
        base = scaf[i]
        if (base[0] != vh) & (base[0] != -1):
            scafPath.append(base[0])
            trans_sc.append([[vh, base[0]]] + [[i, base[1]]])
        if (base[2] != vh) & (base[2] != -1):
            scafPath.append(base[2])
            trans_sc.append([[vh, base[2]]] + [[i, base[3]]])
        if (base[0] == vh) & (base[2] == vh):
            scafPath.append(vh)
        if (base[0] == -1) & (base[2] == -1):
            scafPath.append('-')
        if (base[0] == -1) & (base[2] == vh):
            scafPath.append('e5')
            end5scaf = [vh, i]
        if (base[2] == -1) & (base[0] == vh):
            scafPath.append('e3')
            end3scaf = [vh, i]

    Path['scaf'][vh] = scafPath

    # staple path

    stap = vhToStap[vh]

    for i in range(len(stap)):
        base = stap[i]
        if (base[0] == -1) & (base[2] == vh):
            stapPath.append('e5')
            ends5.append([vh, i])
        if (base[2] == -1) & (base[0] == vh):
            stapPath.append('e3')
            ends3.append([vh, i])
        if (base[0] != vh) & (base[0] != -1):
            stapPath.append(base[0])
            trans.append([[vh, base[0]]] + [[i, base[1]]])
        if (base[2] != vh) & (base[2] != -1):
            stapPath.append(base[2])
            trans.append([[vh, base[2]]] + [[i, base[3]]])
        if (base[0] == vh) & (base[2] == vh):
            stapPath.append(vh)
        if (base[0] == -1) & (base[2] == -1):
            stapPath.append('-')

#    i = 0
#    k = 0
#    while i == 0 and k < len(stapPath):
#        if stapPath[k] != '-':
#            i = 1
#        k += 1
#    if i == 1:
    Path['stap'][vh] = stapPath

# join scaffold and staple
#    print vh, Path['stap'][vh]
#    print vh, Path['scaf'][vh]
    joinedPath = []
    for i in range(len(scafPath)):
        if (scafPath[i] == '-') & (stapPath[i] == '-'):
            joinedPath.append('-')
        elif (scafPath[i] != '-') & (stapPath[i] != '-'):
            joinedPath.append('d')
        elif (scafPath[i] != '-') & (stapPath[i] == '-'):
            joinedPath.append('s')
        else:
            joinedPath.append('o')
    allPath[vh] = joinedPath
#    print vh, allPath[vh]

# ------------ lattice ------------

if args.lattice in ['h', 'hex', 'hexagonal']:
    SQ = False
elif args.lattice in ['s', 'sq', 'square']:
    SQ = True
else:
    raise Exception('Choose lattice type: hexagonal or square')

# ------------- output files ----------

pdb = open(args.output, 'w')
add = open(args.top, 'w')

# ------------- sequence analysis/generation -------------

if args.seq and args.oligs:
    # Get sequence
    seq = ''
    with open(args.seq, 'r') as f:
        for line in f:
            line = line.strip()
            seq += line
    print seq
#    sqdone = []
    # Oligs ends information --> allPath
    with open(args.oligs, 'r') as f:
        for line in f:
            if line[0] != 'S':
                line = line.split(',')
                line[0] = line[0].split('[')
                line[1] = line[1].split('[')
                if line[2][1] != '?':
                    allPath[int(line[0][0])][int(line[0][1][:-1])] = line[2][0]
    #                oligs[int(line[0][0])][int(line[0][1][:-1])] = line[2][0]
                if line[2][-1] != '?':
                    allPath[int(line[1][0])][int(line[1][1][:-1])] = line[2][-1]
    seq_end = None
    compl = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
             'a': 'T', 't': 'A', 'g': 'C', 'c': 'G'}
    [vh, base] = end5scaf
    i = 0
    last1 = False
    last2 = False
    while not seq_end:
        seq_end = checkseq(vh, base)
        if [vh, base] == end3scaf and not seq_end:
            raise Exception('Sequence not accepted')
            seq_end = '-'
        else:
            [vh, base] = nextbase(vh, base, False)
        i += 1
    i = len(seq) - i
    add.write('seq: ' + seq[i:] + seq[:i] + '\n')
elif args.seq or args.oligs:
    raise Exception('Both sequence and oligs are needed')
else:
    add.write('no seq\n')

#ch = (Rows[0] + Cols[0]) % 2

# ---------------- end5 definition (if circular)  ---------------

endtmp = False # indicator whether scaffold is circular
if not end5scaf:
    dl = 10     # length of duplex needed
    endtmp = True
    k = 0
    end = 0
    while not end and k in allPath:
        if (Rows[k] + Cols[k]) % 2:
            c = len(allPath[k]) / 2
            while not end and c < len(allPath[k]) - (dl - 1):
                if allPath[k][c:c + dl] == ['d'] * dl:
                    end = c + dl / 2
                else:
                    c += 1
            c = len(allPath[k]) / 2 - 1
            while not end and c >= 0:
                if allPath[k][c:c + dl] == ['d'] * dl:
                    end = c + dl / 2
                else:
                    c -= 1
            if not end:
                k += 1
        else:
            k += 1
    if not end:
        raise Exception("10bp duplex not found, define scaffold 5'-end mannually")
    end5scaf = [k, end - 1]
    end3scaf = [k, end]

#end5scaf = [0, 45]
#end3scaf = [0, 44]
#[end5scaf, end3scaf] = [end3scaf, end5scaf]
Path['scaf'][end5scaf[0]][end5scaf[1]] = 'e5'
Path['scaf'][end3scaf[0]][end3scaf[1]] = 'e3'

# ------------- particles definition ---------

count = {'D1': 1, 'D2': 2, 'D3': 3, 'D4': 4, 'D5': 5,
         'D6': 6, 'D7': 7, 'T': 7, 'PT': 7, 'S': 1,
         'TT': 7, 'D1T': 1, 'D2T': 2, 'D3T': 3, 'D4T': 4,
         'D5T': 5, 'D6T': 6, 'D7T': 7, 'O': 1, 'OT': 1, 'N': 1}
term = ['T', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'TT', 'D1T', 'D2T',
        'D3T', 'D4T', 'D5T', 'D6T']
if SQ:
    term += ['D7', 'D7T']
    count['PT'] = 8
    count['T'] = 8
    count['D8'] = 8
    count['TT'] = 8

# ------------ scaffold scheme (path) generation -------------

scheme = []
#end = end5scaf
s = 'chain'
k = revers(end5scaf[0], 0)
hm = pathandtype(end5scaf, 'chain')[0]
end = hm[1]
scheme.append(hm[0])
while end[0] != 'end':
    k = revers(end[0], 0)
    [hm, s] = pathandtype(end, s)
    end = hm[1]
    scheme.append(hm[0])

# ---------------- find D7/D8, begin and end for scaffold -------------

for part in scheme:
    k = revers(part[1], 0)
    if part[0] == 'd':      # duplexes
        beg = part[2][0] - part[2][0] % 7 + 7 * (k + 1) / 2 * bool(part[2][0] % 7)
        end = part[2][1] - part[2][1] % 7 - 7 * (k - 1) / 2 * bool(part[2][1] % 7)
        if (end - beg) * k >= 0:    # frf
            part[2] = [(part[2][0], (beg - part[2][0]) * k), range(beg, end + 7 * k, 7 * k),
                       (part[2][1], (part[2][1] - end) * k)]
            if not part[2][0][1]:
                part[2][0] = None
            if not part[2][2][1]:
                part[2][2] = None
        else :      # if duplex is small and without D7/D8 - N particles
            part[2] = range(part[2][0], part[2][1] + k, k)
            part[0] = 'n'
    elif part[0] == 's':    # single-stranded
        part[2] = range(part[2][0], part[2][1] + k, k)


# ------------ print scaffold particles -----------
# ------------ print staple crossovers ------------

if SQ:
    c = 8
    pc = 'D8'
else:
    c = 7
    pc = 'D7'

a_sc = 1    # number of particle printed
nl = 0
atomsh = [None]
anames = [None]
### only 26+10 chains!

l = 1

nls = string.ascii_uppercase + string.digits    # for additional chains
connects = []
needed = []
needst = []
connst = []
scconn = []
scneed = []

sqlines = {}

ssnear = []
ssoligs = []
ssind = []
sstrans1 = []
sstrans2 = []
for v in allPath:
    vh = allPath[v]
#    rev = not((Rows[v] + Cols[v]) % 2)
#    print v, rev
    k = revers(v, 0)
    d = 0
    for i, base in enumerate(vh):
        if base == 'o':
            if Path['stap'][v][i] in ['e3', 'e5']:
                ssoligs.append((v, i))
                ssnear.append(None)
                if (d and k == 1) or (not d and k == -1):
                    ssind.append('e')
                else:
                    ssind.append('b')
                d = not d
            elif Path['stap'][v][i] == v:
                w = False
                if allPath[v][i - k] == 'd':
                    ssoligs.append((v, i))
                    ssnear.append((v, i - k))
                    w = True
                elif allPath[v][i + k] == 'd':
                    ssoligs.append((v, i))
                    ssnear.append((v, i + k))
                    w = True
                if w:
                    if (d and k == 1) or (not d and k == -1):
                        ssind.append('e')
                    else:
                        ssind.append('b')
                    d = not d
            else:
                v1 = Path['stap'][v][i]
                for pair in trans:
                    if pair[0] in [[v, v1], [v1, v]] and i == pair[1][0]:
                        i1 = pair[1][1]
                        if allPath[v1][i1] == 'd':
                            ssoligs.append((v, i))
                            ssnear.append((v1, i1))
                            if (d and k == 1) or (not d and k == -1):
                                ssind.append('e')
                            else:
                                ssind.append('b')
                            d = not d
                        else:
                            sstrans1.append((v, i))
                            sstrans2.append((v1, i1))  ### check v or smth


cnct = open(args.cnct, 'w')
cnct.write('[ distance_restraints ]\n')
cnct.write('; staple crossovers\n')
#cnct.write('; Scaffold-ssolig bonds\n')
for a, part in enumerate(scheme):
    if part[0] == 'd':
        if a == 0:
            if part[2][0]:
                n = 'D' + str(part[2][0][1]) + 'T'
            m = 'TT'
        else:
            if part[2][0]:
                n = 'D' + str(part[2][0][1])
            m = 'T'
        if part[2][0]:
            createatom(part[2][0][0], part[1], n, part[2][0][1], 'b')
            createatom(part[2][1][0], part[1], 'PT', c, 0)
        else:
            createatom(part[2][1][0], part[1], m, c, 'b')
        for i, atom in enumerate(part[2][1][1:-1]):
            createatom(atom, part[1], pc, c, 0)
        if part[2][2]:
            if len(part[2][1]) > 1:
                createatom(part[2][1][-1], part[1], 'PT', part[2][2][1], 0)
            createatom(part[2][2][0], part[1], 'D' + str(part[2][2][1]), 1, 'e')
        else:
            createatom(part[2][1][-1], part[1], 'T', 1, 'e')
    elif part[0] == 's':
            for atom in part[2]:
                createatom(atom, part[1], 'S', 1, 's')
    elif part[0] == 'n':
            if len(part[2]) <= 2:
                for atom in part[2]:
                    createatom(atom, part[1], 'T', 1, 0)
            else:
                createatom(part[2][0], part[1], 'T', 1, 0)
                for atom in part[2][1:-1]:
                    createatom(atom, part[1], 'N', 1, 0)
                createatom(part[2][-1], part[1], 'T', 1, 0)
    else:
        raise Exception("unknown chain type")


# singlestranded 5' !!!
#  check duplexes < 7 bp!

# duplex analysys
# fail if crossover is somewhere else - we doooon't see it


#    print tname, a_sc
# single-stranded staples search

# if crossovers at sspart - FAIL

# oligs length 1 :(
# test crossovers olig-olig


###### !!!! print sstrans!!! if we need

#tn = open(args.tn, 'w')
done = []
while len(done) < len(ssoligs):
    pdb.write('TER\n')
    nl += 1
    i = 0
    while ssoligs[i] in done or ssind[i] != 'b':
        i += 1
    now = ssoligs[i]
    done.append(now)
    [vh, base] = ssoligs[i]
    k = revers(vh, 0)
    createatom(base, vh, 'OT', 1, 's')
    n = ssoligs.index((vh, base))
    if ssnear[n]:
        cnct.write(str(ssnear[n]) + '\t' + str(a_sc - 1) + '\t1\t' + str(l) + '\t1\t0.33\t0.35\t0.37\t0.2\t; ssnear\n')
        l += 1
        add.write('cross ; ' + ' '.join(map(str, [ssnear[n], 0])) + ' ; ' +
              ' '.join(map(str, [a_sc - 1, 0])) + '\n')
    base += k
    while (vh, base) not in ssoligs:
        k = revers(vh, 0)
        createatom(base, vh, 'O', 1, 's')
        if (vh, base) in sstrans1 and (vh, base) not in done:
            n = sstrans1.index((vh, base))
            (vh, base) = sstrans2[n]
            done.append((vh, base))
        elif (vh, base) in sstrans2 and (vh, base) not in done:
            n = sstrans2.index((vh, base))
            (vh, base) = sstrans1[n]
            done.append((vh, base))
        else:
            base += k
    done.append((vh, base))
    createatom(base, vh, 'O', 1, 's')
    n = ssoligs.index((vh, base))
    if ssnear[n]:
        cnct.write(str(ssnear[n]) + '\t' + str(a_sc - 1) + '\t1\t' + str(l) + '\t1\t0.33\t0.35\t0.37\t0.2\n')
        l += 1
        add.write('cross ; ' + ' '.join(map(str, [ssnear[n], 0])) + ' ; ' +
                  ' '.join(map(str, [a_sc - 1, 0])) + '\n')

pdb.close()
#tn.close()


#for i in range(len(connst)):
#    add.write('cross ; ' + ' '.join(map(str, connst[i])) + ' ; ' +
#              ' '.join(map(str, needst[i])) + '\n')
cnct.write('; scaffold crossovers\n')
if endtmp: cnct.write('1\t' + str(a_sc - 1) + '\t1\t' + str(l) + '\t1\t0.33\t0.35\t0.37\t2.5\t; 5\' - 3\' ends\n')
l += 1
for i in range(len(scconn)):
#    write scaffold dist.restr!
    cnct.write(str(scconn[i]) + '\t' + str(scneed[i]) + '\t1\t'
               + str(l) + '\t1\t0.33\t0.35\t0.37\t2.5\n')
    l += 1
    add.write('scaf ; ' + str(scconn[i]) + ' ; ' + str(scneed[i]) + '\n')
cnct.write('\n')
# Print connects file

add.close()

cnct.write('; lattice restraints\n')
if len(connects) != len(needed):
    print 'Smth with connects is wrong'
#for i in range(len(connects)):
#    if connects[i] != connects[i - 1]:
#        cnct.write(str(connects[i]) + '\t' + str(needed[i]) + '\t1\t' + str(l)
#                   + '\t1\t' + '2.0\t2.2\t2.4\t0.5\n')
#        l += 1
cnct.write('\n')

# sequence


# hex and square search

# ---------- lattice restraints ---------

if SQ:
    name = 'D8'
else:
    name = 'D7'
'''
def sqnextline(c, atom):
    global dist
    pp = True
    for k in dist[atom]:
        for ci in c:
            if xy[atom][0][ci] != xy[k][0][ci]:
                pp = False
    if pp:
        n = dist[atom].index(k)
        del dist[atom][n]
        n = dist[k].index(atom)
        del dist[k][n]
        return k

def addsq(atom, atom2, key):
    f = sqlines.get(key, [])
    res = []
    if atom not in f:
        res.append(atom)
    if atom2 not in f:
        res.append(atom2)
    return res


def twsearch(c1, c2):
    if c1[0] == c2[0]:
        if c1[1] == c2[1]:
            return (c1[0], c1[1], None)
        elif c1[2] == c2[2]:
            return (c1[0], None, c1[2])
    elif c1[1:] == c2[1:]:
        return (None, c1[1], c1[2])
    return False


def longsearch(a, b):
    a = int(a)
    b = int(b)
    (vh, i) = atomsh[a]
    (z, i1) = atomsh[b]
    c = -1
    n1 = None
    n2 = None
    while Path['scaf'][vh][i + c] != '-'  and Path['scaf'][z][i + c] != '-' and not n1:
        if Path['stap'][vh][i + c] == z:
            return True
        c -= 1
    c = 1
    while Path['scaf'][vh][i + c] != '-' and Path['scaf'][z][i + c] != '-' and not n2: ### != vh, z
        if Path['stap'][vh][i + c] == z:
            return True
        c += 1
#    return False
    return True


def asc2coords(atom):
    atom = int(atom)
    return [atomsh[atom][1], Rows[atomsh[atom][0]], Cols[atomsh[atom][0]]]

pymol.finish_launching()
pymol.cmd.load(args.output)
pymol.stored.list = []
pymol.cmd.iterate("(resn " + name + "+T+TT+PT)", "stored.list.append(resi)")
atoms = pymol.stored.list
dist = {}
for atom in atoms:
    pymol.stored.list = []
    pymol.cmd.iterate("(resn " + name + "+T+TT+PT and not resi " + str(atom) + " within 23.1 of resi " + str(atom) + ")", "stored.list.append(resi)")
    if pymol.stored.list:
        dist[atom] = pymol.stored.list

def checkthird(center, new, d):
    if new[0] == center[0]:
        if d == 1: c1 = 2
        elif d == 2: c1 = 1
        if center[d] == new[d]:
            return True
        elif center[c1] == new[c1] and not SQ:
            return True
    return False

if SQ:
    [a, b, c] = ['4.54', '4.6', '4.63']
else:
    [a, b, c] = ['3.68', '3.75', '3.77']

cnct.write('; lattice\n')

#print dist['80']

if True:
    tmp = dist.keys()
    done = []
    while tmp:
        sdone = []
        first = tmp[0]
        done.append(first)
        fcoords = asc2coords(first)
        for second in dist[first]:
#            if first == '83': print second, '2'
            sdone.append(second)
            scoords = asc2coords(second)
            if scoords[0] == fcoords[0]:
                if scoords[1] == fcoords[1]:
                    sim = 1
                elif scoords[2] == fcoords[2]:
                    sim = 2
                else:
                    print '?'
            for third in dist[first]:
                if third not in done and third not in sdone:
                    tcoords = asc2coords(third)
                    if checkthird(fcoords, tcoords, sim):
                        if longsearch(first, third) and longsearch(first, second):
                            cnct.write(third + '\t' + second + '\t1\t' + str(l) + '\t1\t' + a + '\t' + b + '\t' + c + '\t1.5\n')
                            l += 1
            for third in dist[second]:
                if third not in done:
#                    if first == '83': print third, '3', third in done
                    tcoords = asc2coords(third)
                    if checkthird(scoords, tcoords, sim):
#                        if first == '83': print third, '3', third in done, 'checked'
                        if longsearch(second, third) and longsearch(first, second):
                            if first == '83': print third, '3', third in done, 'and where?'
                            cnct.write(third + '\t' + first + '\t1\t' + str(l) + '\t1\t' + a + '\t' + b + '\t' + c + '\t1.5\n')
                            l += 1
        for second in dist[first]:
            n = dist[second].index(first)
            del dist[second][n]
        tmp.pop(0)

pymol.cmd.quit()
'''

def crosscheck(one, two, thr):
    '''Check if crossovers in the triple'''
    c = -1
    ok = True
    if SQ: lim = 29
    else:   lim = 100
    for pair in [[one, two], [two, thr]]:
        ok1 = False
        ok2 = False
        [[vh1, b1], [vh2, b2]] = pair
        while Path['scaf'][vh1][b1 + c] != '-' and Path['scaf'][vh2][b2 + c] != '-' and c > lim * -1:
#            print Path['stap'][vh1][b1 + c]
            if Path['stap'][vh1][b1 + c] == vh2:
                ok1 = True
            c -= 1
        c = 1
        while Path['scaf'][vh1][b1 + c] != '-' and Path['scaf'][vh2][b2 + c] != '-' and c < lim:
#            print Path['stap'][vh1][b1 + c]
            if Path['stap'][vh1][b1 + c] == vh2:
                ok2 = True
            c += 1
        ok = ok & (ok1 | ok2)
    return ok

def findtriples(base, coords):
    (row, col) = base
    res = []
    cpair = [None, None]
    rpair1 = [None, None]
    rpair2 = [None, None]
    rpair3 = [None, None]
    rpair4 = [None, None]
    if SQ:
        for i, atom in enumerate(coords):
            if atom == (row, col + 1): cpair[0] = i
            elif atom == (row, col + 2): cpair[1] = i
            elif atom == (row + 1, col): rpair1[0] = i
            elif atom == (row + 2, col): rpair1[1] = i
    else:
        for i, atom in enumerate(coords):
            if atom == (row, col + 1): cpair[0] = i
            elif atom == (row, col + 2): cpair[1]= i
            if (row + col) % 2:
                if atom == (row + 1, col): rpair3[0] = i; rpair4[0] = i
                elif atom == (row + 1, col + 1): rpair3[1] = i
                elif atom == (row + 1, col - 1): rpair4[1] = i
            else:
                if atom == (row, col + 1): rpair1[0] = i
                elif atom == (row + 1, col + 1): rpair1[1] = i
                elif atom == (row, col - 1): rpair2[0] = i
                elif atom == (row + 1, col - 1): rpair2[1] = i
    for pair in [cpair, rpair1, rpair2, rpair3, rpair4]:
        if None not in pair:
            res += [pair]
    return res

if SQ:
    [a, b, c] = ['4.54', '4.6', '4.63']
else:
    [a, b, c] = ['3.68', '3.75', '3.77']

triples = {}
trcoords = {}
for i, atom in enumerate(atomsh):
    if atom and anames[i] in [name, 'T', 'TT', 'PT']:
        x = atom[1]
        trcoords[x] = trcoords.get(x, []) + [(Rows[atom[0]], Cols[atom[0]])]
        triples[x] = triples.get(x, []) + [i]
for x in trcoords:
    vhs = trcoords[x]
    for i, base in enumerate(vhs):
        tr = findtriples(base, vhs)
        for pair in tr:
            one = atomsh[triples[x][i]]
            two = atomsh[triples[x][pair[0]]]
            thr = atomsh[triples[x][pair[1]]]
            if crosscheck(one, two, thr):
                cnct.write(str(triples[x][i]) + '\t' + str(triples[x][pair[1]]) + '\t1\t' + str(l) + '\t1\t' + a + '\t' + b + '\t' + c + '\t1.5\n')
                l += 1

cnct.close()
