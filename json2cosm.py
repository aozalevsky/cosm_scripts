#!/usr/bin/env python
# encoding: utf-8

import argparse
import string
import copy
import math
import util

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True,
                    action='store', dest='input',
                    help='input cadnano json file')
parser.add_argument('-o', '--output', required=True,
                    action='store', dest='output',
                    help='output pdb file')
parser.add_argument('-r', '--restr', required=True,
                    action='store', dest='cnct',
                    help='output distance/angle restraints')
parser.add_argument('-t', '--top', required=True,
                    action='store', dest='top',
                    help='output file for ccg --> cg / fullatom convertation')
parser.add_argument('-l', '--lattice', required=True,
                    action='store', dest='lattice',
                    help='lattice type (hexagonal / square)')
parser.add_argument('-m', '--map', required=True,
                    action='store', dest='map',
                    help='2d coords map file')
parser.add_argument('--seq',
                    action='store', dest='seq',
                    help='scaffold sequence')
parser.add_argument('--oligs',
                    action='store', dest='oligs',
                    help='staple sequences exported csv')
args = parser.parse_args()


# do not 3-5 restr if additional oligs
# lattice crossovers - not only near
# other constants for olig-ssnear (hexcl)


HELIXDIST = 22.0


def honeycomb(row, col, dist=HELIXDIST):
    """Cadnano honeycomb lattice into cartesian"""
    # k = distance between helixes
    ANGLE = 30.0
    x = col * dist * math.cos(math.radians(ANGLE))
    if (row % 2) != 0:
        y = (row * 3 + 1 - col % 2) * dist * math.sin(math.radians(ANGLE))
    else:
        y = (row * 3 + col % 2) * dist * math.sin(math.radians(ANGLE))
    return x, y


def square(row, col, dist=HELIXDIST):
    """Cadnano square lattice into cartesian"""
    x = col * dist
    y = row * dist
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
    res = None
#    print Path['scaf'][3][5]
#    print cross_base, cross_vh
    cross_base = end[1] + k * cross_base + k  # + (k + 1)/2
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
    # works two times in sequence search:
    # 1st iteration for next sequence negin point
    # 2nd iterantion for sequence checking
    global last1, last2
    if seq:
        lastc = last2
    else:
        lastc = last1
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
    # k = revers(vh, 0) # Commented during code cleaning
    i = 0
    while [vh, base] != beg or not i:  # and Path['scaf'][vh][base] != 'e3':
        if endtmp and Path['scaf'][vh][base] == 'e3':
            return True
        if seqPath[vh][base] != '-':  # in ['A', 'T', 'G', 'C']:
            if seq[i].upper() != compl[seqPath[vh][base]]:
                return False
        [vh, base] = nextbase(vh, base, True)
        i += 1
    return True


def findnextnot(end, pat):
    global k
#    print end, pat
    if k == 1:
        path = allPath[end[0]][end[1] + 1:]
    else:
        path = allPath[end[0]][:end[1]]
        path.reverse()
    if path:
        path = [pat] + path
        i = 0
        cross_base = None
        while i < len(path) - 1 and not cross_base:
            i += 1
            if path[i] not in [pat, '-']:
                cross_base = i

    #   cross_vh = path[cross_base]
        if cross_base:
            #   cross_base += 1
            cross_base = end[1] + k * cross_base
#            print cross_base
#            if k == -1:
#                cross_base -= 1
            return [end, [end[0], cross_base]]
        else:
            return None
    else:
        return None


def outscheme(a, b):
    mapf.write('L ' + str(a) + ':' + str(vhNums.index(b)) + '\n')


def pathandtype(end, s):
    """Get next duplex end or vh (chain) change"""
    global k
    if s == 'change':
        new_chain = get_scaffold_path([end[0], end[1] - k])
    else:
        new_chain = get_scaffold_path(end)
    new_change = findnextnot(end, allPath[end[0]][end[1]])
    if new_chain[1] is None:
        new_end = new_change[1]
        s = 'change'
        l = new_change[1][1] - k
        return [([allPath[end[0]][end[1]], end[0], [end[1], l]], new_end), s]
    if new_change:
        if new_chain[0][1] * k <= (new_change[1][1] - k) * k:
            new_end = new_chain[1]
            s = 'chain'
            l = new_chain[0][1]
            outscheme(new_chain[0][1], new_chain[0][0])
            if new_chain[1][0] != 'end':
                outscheme(new_chain[1][1], new_chain[1][0])

        else:
            new_end = new_change[1]
            s = 'change'
            l = new_change[1][1] - k
    else:
        new_end = new_chain[1]
        s = 'chain'
        l = new_chain[0][1]
        outscheme(new_chain[0][1], new_chain[0][0])
        if new_chain[1][0] != 'end':
            outscheme(new_chain[1][1], new_chain[1][0])
    return [([allPath[end[0]][end[1]], end[0], [end[1], l]], new_end), s]


def prconnect(base1, base2):
    """Print crossover restraints"""
    global l
    template = "{0[0]:<6s}{0[1]:>5d}{0[2]:>5d}"
    t = tuple(['CONECT', base1, base2])
    outpdbc.append(template.format(t) + '\n')
    if not (isinstance(base2, int) and isinstance(base1, int)):
        raise Exception('ADMIN: PY1. Error with staple restraints')


def findconnects(vh, j, number, length, end):
    """Find crossover restraints"""
    global connects, needed, add_tmp
    if end == 's':
        return None
    # [x, y] = [0, 0]
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
                    add_tmp.append([number, 0, connst[n][0], connst[n][1]])
                elif i % HC == HC - 1:
                    needst[n] = (number, dif)
                    if k == -1:
                        add_tmp.append([number, 1, connst[n][0], connst[n][1]])
                    else:
                        add_tmp.append(
                            [number, length - 1, connst[n][0], connst[n][1]])
                else:
                    needst[n] = (number + dif, 0)
                    add_tmp.append(
                        [number + dif, 0, connst[n][0], connst[n][1]])
            if [vh, i] in needed:
                n = needed.index([vh, i])
                if i == j:
                    needed[n] = number
                elif i % HC == HC - 1:
                    if k == 1:
                        needed[n] = number + 1
                    else:
                        needed[n] = number
                else:
                    ind.append(dif)
                    needed[n] = number + dif
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
                        elif i % HC == HC - 1:
                            if k == -1:
                                connst.append((number, 1))
                            else:
                                connst.append((number, length - 1))
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
                            if dif not in ind:
                                ind.append(dif)
                            connects.append(number + dif)
                            needed.append([pair[0][1], pair[1][1]])

    if ind:
        return ind
    return None


def checkreg(pair):
    [[vh1, b1], [vh2, b2]] = pair
    if not (isinstance(b1, int) and isinstance(b2, int)):
        return True
    c = -1
    ok1 = False
    ok2 = False
    end = ['-']
    if endtmp:
        end += ['e5', 'e3']
    while (Path['scaf'][vh1][b1 + c] not in end and
            Path['scaf'][vh2][b2 + c] not in end and
            c > LIM * -1):
        if Path['stap'][vh1][b1 + c] == vh2:
            ok1 = True
        c -= 1
    c = 1

    while (Path['scaf'][vh1][b1 + c] not in end and
            Path['scaf'][vh2][b2 + c] not in end and
            c < LIM):
        if Path['stap'][vh1][b1 + c] == vh2:
            ok2 = True
        c += 1
    return ok1 | ok2


def crosscheck(one, two, thr):
    '''Check if crossovers in the range'''
    ok = True
    # print one, two, thr
    for pair in [[one, two], [two, thr]]:
        ok = ok & checkreg(pair)
    return ok


def TtoB(name, atom1):
    if atom1 == a_sc - 1 or atom1 == 1:
        return name
    atom2 = atom1 + 1
    atom1 = atomsh[atom1]
    atom2 = atomsh[atom2]
    if checkreg([atom1, atom2]):
        return name
    else:
        if name == 'T':
            return 'B'
        else:
            return 'B' + name[1:]


def check_circular(vh, i):
    global valStap
    k = revers(vh, 1)
    end5 = (vh, i)
    end3 = False
    nextbase = None
    prevcr = False
    while not end3 and nextbase != end5:
        valStap[vh][i] = '-'
        if nextbase:
            (vh, i) = nextbase
        else:
            (vh, i) = end5
        k = revers(vh, 1)
        nextbase = None
#        print vh, i, Path['stap'][vh][i]
        if Path['stap'][vh][i] in [vh, 'e5'] or prevcr:
            nextbase = (vh, i + k)
            prevcr = False
        else:
            if not prevcr:
                for pair in trans:
                    if pair[0][0] == vh and pair[1][0] == i:
                        nextbase = (pair[0][1], pair[1][1])
                        prevcr = True
            if not nextbase:
                if Path['stap'][vh][i] == 'e3':
                    end3 = True
                    valStap[vh][i] = '-'
                else:
                    raise Exception(
                        'ADMIN: PY1. Smth wrong with staple crossover paths')
    if nextbase == end5:
        raise Exception('ADMIN: PY1. Smth wrong with staple crossover paths')


def staplesends(vh, i, number, length):
    """Find staple ends through path"""
    global needst, connst, add_ends
    rev = (Rows[vh] + Cols[vh]) % 2
    if not isinstance(i, int):
        return False
    if not rev:
        k = 1
    else:
        k = -1
    l = 0
    while l < length:
        a = Path['stap'][vh][i]
        if a == 'e5':
            add_ends.append(('e5', vh, i))
            check_circular(vh, i)
        elif a == 'e3':
            add_ends.append(('e3', vh, i))
        i += k
        l += 1


def scaffoldcross(vh, i, number):
    """Find and print scaffold crossovers"""
    global scconn, scneed  # check for inserts!
    rev = (Rows[vh] + Cols[vh]) % 2
    if not rev:
        k = 1
    else:
        k = -1
    if (vh, i) in scneed:  # test this
        n = scneed.index((vh, i))
#        t = atomsh[scconn][n]
#        if number
        scneed[n] = number
    elif (vh, i) in scconn:
        n = scconn.index((vh, i))
        scconn[n] = number
    else:
        n = None
        for pair in trans_sc:
            if pair[0][0] == vh and pair[1][0] == i:
                if Path['scaf'][vh][i + k] not in [vh, '-']:
                    if Path['scaf'][vh][i - k] in [vh, '-']:
                        if (vh, i + k) in atomsh:
                            n = atomsh.index((vh, i + k))
                        if n != number - 1:
                            scneed.append((vh, i + k))
                            scconn.append(number)
                    else:
                        raise Exception(
                            'ADMIN: PY1. Error with scaffold '
                            'restraints (2bp duplex)')
                elif Path['scaf'][vh][i - k] not in [vh, '-']:
                    if (vh, i - k) in atomsh:
                        n = atomsh.index((vh, i - k))
                    if n != number - 1:
                        scconn.append((vh, i - k))
                        scneed.append(number)


def outatom(atom):
    """Print atom in pdb"""
    global a_sc
    if len(atom) > 1:
        template = ("{0[0]:<6s}{0[1]:>5d}  {0[2]:<4s}{0[3]:<3s} "
                    "{0[4]:>1s}{0[5]:>4d}    "
                    "{0[6]:>8.3f}{0[7]:>8.3f}{0[8]:>8.3f}"
                    "{0[9]:>6.2f}{0[10]:>6.2f}            {0[11]:<2s}"
                    )
        t = tuple(atom[:-2])
        pdb.write(template.format(t) + '\n')
    else:
        pdb.write(atom[0])


def addtoout(aname, tname, nlname, ox, oy, z, chain, length, mod):
    global a_sc, outpdb, atomsh, anames, outmap
    if tname:
        t1 = [aname] + [tname] + [nlname]
        t2 = [ox, oy] + [z * 3.4] + [1.00, 0.00] + [aname, chain, z]
        if mod == 'I':
            outpdb.append(['ATOM', a_sc] + t1 + [a_sc] + t2)
            outmap.append(['T', tname, a_sc, chain, z])
            staplesends(chain, z, a_sc, length)
            scaffoldcross(chain, z, a_sc)
            atomsh.append((chain, z))
            anames.append(tname)
            a_sc += 1
            t1[1] = 'N'
            outpdb.append(['ATOM', a_sc] + t1 + [a_sc] + t2)
            outmap.append(['I', 'N', a_sc, chain, z])

            atomsh.append((chain, z))
            anames.append('N')
            a_sc += 1
        elif mod == 'Ie':
            savename = tname
            outpdb.append(['ATOM', a_sc] + t1 + [a_sc] + t2)
            outmap.append(['I', tname, a_sc, chain, z])
            atomsh.append((chain, z))
            anames.append('N')
            a_sc += 1
            t1[1] = savename
            outpdb.append(['ATOM', a_sc] + t1 + [a_sc] + t2)
            outmap.append(['T', tname, a_sc, chain, z])
            staplesends(chain, z, a_sc, length)
            scaffoldcross(chain, z, a_sc)
            atomsh.append((chain, z))
            anames.append(savename)
            a_sc += 1
        elif not mod:
            outpdb.append(['ATOM', a_sc] + t1 + [a_sc] + t2)
            outmap.append(['T', tname, a_sc, chain, z])
            staplesends(chain, z, a_sc, length)
            scaffoldcross(chain, z, a_sc)
            atomsh.append((chain, z))
            anames.append(tname)
            a_sc += 1
        elif mod == 'D':
            allPath[chain][z] = 'DELE'
        else:
            raise Exception('ADMIN: PY1. Wrong modification type')


def createatom(z, chain, tname, length, end, modif):
    """Print particles in the interval -- calculate coords etc"""
    global a_sc, dif
    [row, col] = [Rows[chain], Cols[chain]]
    if not SQ:
        [ox, oy] = list(honeycomb(row, col))
    else:
        [ox, oy] = list(square(row, col))
    if a_sc == 1:
        dif = None
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
    if set(inserts) and set(delet):
        raise Exception(
            'USER: Insertion and deletion in the same place.')
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
    savedif = False
    if tname == 'N':
        if delet:
            addtoout(aname, 'N', nlname, ox, oy, z, chain, 1, 'D')
        if inserts:
            addtoout(aname, 'N', nlname, ox, oy, z, chain, 1, 'I')
    else:
        if dif and tname not in term:
            tname = 'PT'
        if dif and tname in term and tname[0] == 'T':
            if not inserts and not delet:
                tname = 'T'
            else:
                length = int(tname[1])
                tname = 'T'
                savedif = True
        dif = findconnects(chain, z, a_sc, length, end)
        if dif:
            if tname not in term:
                tname = 'PT'
            if tname in term and tname[0] == 'T':
                tname = 'T'
            if z in inserts:
                M = 'I'
            elif z in delet:
                M = 'D'
            else:
                M = None
            addtoout(aname, tname, nlname, ox, oy, z, chain, 1, M)
            for i in range(1, length):
                if z + i * k in inserts:
                    M = 'I'
                elif z + i * k in delet:
                    M = 'D'
                else:
                    M = None
                addtoout(aname, 'N', nlname, ox, oy, z + i * k, chain, 1, M)
        elif inserts or delet:
            if not instype:
                if tname not in term:
                    tname = 'PT'
                if tname in term and tname[0] == 'T':
                    tname = 'T'
                if z in inserts:
                    M = 'I'
                elif z in delet:
                    M = 'D'
                else:
                    M = None
                if length > 1:
                    addtoout(aname, tname, nlname, ox, oy, z, chain, length, M)
                elif M:
                    addtoout(aname, tname, nlname, ox, oy,
                             z, chain, length, M + 'e')  # first insert
                else:
                    addtoout(
                        aname, tname, nlname, ox, oy, z, chain, length, None)
                for i in range(1, length):
                    if z + i * k in inserts:
                        M = 'I'
                    elif z + i * k in delet:
                        M = 'D'
                    else:
                        M = None
                    addtoout(
                        aname, 'N', nlname, ox, oy, z + i * k, chain, 1, M)
                dif = True
            else:
                if not savedif:
                    length = count[tname]
                if z in delet and tname == 'T1':
                    raise Exception(
                        'USER: Deletion at terminal base pair')  # because of T1
                    MT = 'D'
                if z in inserts:
                    MT = 'I'
                else:
                    MT = None
                if instype == 't':
                    addtoout(aname, 'T', nlname, ox, oy, z, chain, length, MT)
                    for i in range(1, length):
                        M = None
                        if z + i * k in inserts:
                            M = 'I'
                        elif z + i * k in delet:
                            M = 'D'
                        addtoout(
                            aname, 'N', nlname, ox, oy, z + i * k, chain, 1, M)
                else:
                    for i in range(1, length):
                        M = None
                        if z + (i - length) * k in inserts:
                            M = 'I'
                        elif z + (i - length) * k in delet:
                            M = 'D'
                        addtoout(aname, 'N', nlname, ox, oy, z + (
                            i - length) * k, chain, 1, M)
                    if MT:
                        addtoout(
                            aname, 'T', nlname, ox, oy, z,
                            chain, length, MT + 'e')
                    else:
                        addtoout(
                            aname, 'T', nlname, ox, oy, z,
                            chain, 1, None)  # !!!
        else:
            addtoout(aname, tname, nlname, ox, oy, z,
                     chain, length, None)  # quite normal atom


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

# ------------- .json analysis ----------------
outpdb = []
outpdbc = []
outmap = []


obj = util.load_json(args.input)
util.check_json_version(obj)

try:
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
    loop = {}
    skip = {}
    end5scaf = None
    end3scaf = None
    RCvh = {}

    for strand in strands:
        num = strand["num"]
        if num in vhNums:
            raise Exception("USER: Double strand names in json file")
        vhNums.append(num)
        Rows[num] = strand["row"]
        Cols[num] = strand["col"]
        scaf = strand["scaf"]
        stap = strand["stap"]
        vhToScaf[num] = scaf
        vhToStap[num] = stap
        loop[num] = strand["loop"]
        skip[num] = strand["skip"]
        RCvh[(strand["row"], strand["col"])] = num

    scaf5 = []
    scaf3 = []
    coords = []

    Path = {'scaf': {}, 'stap': {}}
    ends5 = []
    ends3 = []
    allPath = {}
    seqPath = {}

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
        seqPath[vh] = len(joinedPath) * ['-']

except Exception as e:
    print(e)
    raise Exception('USER: Error in input json file')

print 'end5scaf', end5scaf
print 'end3scaf', end3scaf

for v in allPath:
    vh = allPath[v]
    if 'o' in vh:
        raise Exception('USER: Single-stranded staple')


# ------------- loop --------------
add = open(args.top, 'w')
insert = {}
deletion = {}
mapf = open(args.map, 'w')

for vh in loop:
    for i, n in enumerate(loop[vh]):
        # if n != 0:  # Useless check
        if n > 1:
            raise Exception(
                'USER: More than one base inserted in one place')
        else:
            insert[(vh, i)] = n

for vh in skip:
    for i, n in enumerate(skip[vh]):
        # if n != 0:
        if abs(n) > 1:
            raise Exception(
                'USER: More than one base deleted in one place')
        elif n != 0:
            deletion[(vh, i)] = abs(n)
            mapf.write('D ' + str(i) + ':' + str(vhNums.index(vh)) + '\n')

# ------------ lattice ------------

if args.lattice in ['h', 'hex', 'hexagonal']:
    SQ = False
elif args.lattice in ['s', 'sq', 'square']:
    SQ = True
else:
    raise Exception('ADMIN: Wrong lattice type')

if SQ:
    LIM = 36 + 7
    HC = 8
else:
    LIM = 21 + 6
    HC = 7

# ------------- output files ----------


# ch = (Rows[0] + Cols[0]) % 2

# ---------------- end5 definition (if circular)  ---------------

endtmp = True  # indicator whether scaffold is circular
if not end5scaf:
    if SQ:
        dl = 36
    else:
        dl = 21     # length of duplex needed
    k = min(vhNums)
    end = 0
    endtmp = False
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
        end5scaf = [k, end - revers(k, 1)]
        end3scaf = [k, end]
    if not end:
        p = min(vhNums)
        while not end and p in allPath:
            b = 1
            while not end and b < len(Path['scaf'][p]):
                if allPath[p][b - 1: b + 2] == ['d'] * 3:
                    end = True
                    end5scaf = [p, b]
                    end3scaf = [p, b + revers(p, 1)]
                b += 1
            p += 1
        if not end:
            raise Exception("ADMIN: PY1. 5'-end search error")
Path['scaf'][end5scaf[0]][end5scaf[1]] = 'e5'
Path['scaf'][end3scaf[0]][end3scaf[1]] = 'e3'


# ------------- sequence analysis/generation -------------

compl = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
         'a': 'T', 't': 'A', 'g': 'C', 'c': 'G'}

if args.seq and args.oligs:
    # Get sequence
    seq = ''
    try:
        with open(args.seq, 'r') as f:
            for line in f:
                line = line.strip()
                seq += line
        last1 = False
        last2 = False
    except:
        raise Exception('USER: Error in sequence input file')
    for letter in seq:
        if letter not in compl:
            raise Exception('USER: Error in sequence input file')

    try:
        with open(args.oligs, 'r') as f:
            for line in f:
                if line[0] != 'S':
                    line = line.split(',')
                    line[0] = line[0].split('[')
                    line[1] = line[1].split('[')
                    if line[2][1] != '?':
                        seqPath[int(line[0][0])][
                            int(line[0][1][:-1])] = line[2][0]
                    if line[2][-1] != '?':
                        seqPath[int(line[1][0])][
                            int(line[1][1][:-1])] = line[2][-1]
    except:
        raise Exception('USER: Error in staple csv file')
    for vh in seqPath:
        for i in seqPath[vh]:
            if i not in compl and i != '-':
                raise Exception('USER: Error in staple csv file')

    if endtmp:
        add.write('seq: ' + seq + '\n')
        if not checkseq(end5scaf[0], end5scaf[1]):
            raise Exception('USER: Sequence does not match staple list')
    else:
        seq_end = None
        [vh, base] = end5scaf
        i = 0

        while [vh, base] != end3scaf:
            seq_end = checkseq(vh, base)
            [vh, base] = nextbase(vh, base, False)
            if seq_end:
                break
            i += 1
        if not seq_end:
            raise Exception('USER: Sequence does not match staple list')
            seq_end = '-'
        i = len(seq) - i
        add.write('seq: ' + seq[i:] + seq[:i] + '\n')
elif args.seq or args.oligs:
    raise Exception('USER: Both scaffold and staple sequences are needed')
else:
    add.write('no seq\n')


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
if SQ:
    count['PT'] = 8
    count['T'] = 8
    count['B'] = 8
    count['BT'] = 8
    count['H'] = 8
    count['TT'] = 8

# ------------ scaffold scheme (path) generation -------------

scheme = []
s = 'chain'
k = revers(end5scaf[0], 0)
mapf.write('L ' + str(end5scaf[1]) + ':' + str(
    vhNums.index(end5scaf[0])) + '\n')
[hm, s] = pathandtype(end5scaf, s)
end = hm[1]
scheme.append(hm[0])
while end[0] != 'end':
    k = revers(end[0], 0)
    [hm, s] = pathandtype(end, s)
    end = hm[1]
    scheme.append(hm[0])


valPath = copy.deepcopy(Path['scaf'])
valStap = copy.deepcopy(Path['stap'])

# -------- more than one scaffolds validation -------------

for p in scheme:
    pb = max(p[2])
    pe = min(p[2])
    for i in range(pe, pb + 1):
        valPath[p[1]][i] = '-'

for vh in valPath:
    for i in valPath[vh]:
        if i != '-':
            raise Exception('USER: More than one scaffold chain')
# ---------------------------------------------------------

if SQ:
    dp = 0
else:
    dp = 0
for i, part in enumerate(scheme):
    scheme[i][2] = [part[2][0] + dp, part[2][1] + dp]
for vh in allPath:
    allPath[vh] = ['-'] * dp + allPath[vh]
    Path['scaf'][vh] = ['-'] * dp + Path['scaf'][vh]
    Path['stap'][vh] = ['-'] * dp + Path['stap'][vh]

# ---------------- find H, begin and end for scaffold -------------

new_scheme = []
for part in scheme:
    k = revers(part[1], 0)
    if part[0] == 'd':      # duplexes
        ins1 = []
        ins2 = []
        ins = []
        dele1 = []
        dele2 = []
        dele = []
        beg = part[2][0] - part[2][0] % HC + \
            HC * (k + 1) / 2 * bool(part[2][0] % HC)
        end = part[2][1] - part[2][1] % HC - \
            HC * (k - 1) / 2 * bool(part[2][1] % HC)
        isH = False
        for base in range(part[2][0], part[2][1] + k, k):
            if base % HC:
                isH = True
        if (end - beg) * k >= 0 and isH:
            part[2] = [(
                part[2][0], (
                    beg - part[2][0]) * k), range(beg, end + HC * k, HC * k),
                (part[2][1], (part[2][1] - end) * k)]

            # find insertions
            for i in range(part[2][0][0], part[2][1][0], k):    # first term
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
                    for l in range(b, b + k * HC, k):
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
                    for l in range(b, b + k * HC, k):
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
                raise Exception(
                    'USER: Insertion or deletion in single-stranded region')

# ------------ print scaffold particles -----------
# ------------ print staple crossovers ------------

pc = 'H'
a_sc = 1    # number of particle printed
nl = 0
atomsh = [None]
anames = [None]
add_tmp = []
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
    k = revers(v, 0)
    d = 0
    for i, base in enumerate(vh):
        if base == 'o':
            raise Exception('USER: Single-stranded staple')
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
                            sstrans2.append((v1, i1))  # check v or smth

cnct = open(args.cnct, 'w')
cnct.write('[ distance_restraints ]\n')
cnct.write('; staple crossovers\n')

add_ends = []

for a, part in enumerate(scheme):
    if part[0] == 'd':
        # first term

        if part[2][0]:
            n = 'T' + str(part[2][0][1])
        m = 'T'
        tmpins = 0
        tmpdel = 0
        if part[2][0]:
            createatom(part[2][0][0], part[1], n, part[2][
                       0][1], 'b', ('t', part[2][0][2], part[2][0][3]))
            if len(part[2][1]) <= 1:
                if part[2][2]:
                    createatom(part[2][1][0], part[
                               1], 'PT', 1, 0, (part[3][0], part[4][0]))
                else:
                    createatom(part[2][1][0], part[
                               1], 'T', 1, 0, (part[3][0], part[4][0]))
            else:
                createatom(part[2][1][0], part[
                           1], 'PT', HC, 0, (part[3][0], part[4][0]))

        else:
            createatom(
                part[2][1][0], part[1], m, HC, 'b', (part[3][0], part[4][0]))

        # center
        for i, atom in enumerate(part[2][1][1:-1]):
            createatom(atom, part[1], pc, HC, 0, (
                part[3][i + 1], part[4][i + 1]))

        # last term

        if part[2][2]:
            if len(part[2][1]) > 1:
                createatom(part[2][1][-1], part[1], 'PT', part[
                           2][2][1], 0, (part[3][-1], part[4][-1]))
            createatom(part[2][2][0], part[1], 'T' + str(
                part[2][2][1]), 1, 'e', ('e', part[2][2][2], part[2][2][3]))
        else:
            if len(part[2][1]) > 1:
                createatom(part[2][1][-1], part[
                           1], 'T', 1, 'e', (part[3][-1], part[4][-1]))
    elif part[0] == 's':
        if a == 0:
            createatom(part[2][0], part[1], 'S', 1, 's', False)
            if len(part[2]) > 1:
                for atom in part[2][1:]:
                    createatom(atom, part[1], 'S', 1, 's', False)
        else:
            for atom in part[2]:
                createatom(atom, part[1], 'S', 1, 's', False)
    elif part[0] == 'n':
        tmpins = False
        tmpdel = False
        if len(part[2]) <= 2:
            for atom in part[2]:
                if atom in part[4]:
                    raise Exception('USER: Deletion at terminal base pair')
                if atom in part[3]:
                    createatom(atom, part[1], 'T', 1, 0, (1, 0))
                else:
                    createatom(atom, part[1], 'T', 1, 0, (0, 0))
        else:
            if part[2][0] in part[4]:
                raise Exception('USER: Deletion at terminal particle')
            if part[2][-1] in part[4]:
                raise Exception('USER: Deletion at terminal particle')
            if part[2][0] in part[3]:
                createatom(part[2][0], part[1], 'T', 1, 0, (1, 0))
            else:
                createatom(part[2][0], part[1], 'T', 1, 0, (0, 0))
            for atom in part[2][1:-1]:
                tmpins = 0
                tmpdel = 0
                if atom in part[3]:
                    tmpins = 1
                if atom in part[4]:
                    tmpdel = 1
                createatom(atom, part[1], 'N', 1, 0, (tmpins, tmpdel))
            if part[2][-1] in part[3]:
                createatom(part[2][-1], part[1], 'T', 1, 0, (1, 0))
            else:
                createatom(part[2][-1], part[1], 'T', 1, 0, (0, 0))
    else:
        raise Exception("ADMIN: PY1. Unknown chain type")


# singlestranded 5' !!!
#  check duplexes < 7 bp!

# single-stranded staples search

# if crossovers at sspart - FAIL

# oligs length 1 :(
# test crossovers olig-olig


# !!!! print sstrans!!! if we need

# tn = open(args.tn, 'w')
# sequence


# hex and square search

# -------- check circular staples --------

for vh in valStap:
    for i in valStap[vh]:
        if i != '-':
            raise Exception('USER: Circular staple')

# ---------------------------------------

done = []
num = len(outpdb) + 1
while len(done) < len(ssoligs):
    outpdb.append(['TER\n'])
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
        cnct.write(str(ssnear[n]) + '\t' + str(a_sc - 1) + '\t1\t' + str(
            l) + '\t1\t0.33\t0.35\t0.37\t0.2\t; ssnear\n')
        l += 1
        add.write('cross ; ' + ' '.join(map(str, [ssnear[n], 0])) + ' ' +
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
        cnct.write(str(ssnear[n]) + '\t' + str(
            a_sc - 1) + '\t1\t' + str(l) + '\t1\t0.33\t0.35\t0.37\t0.2\n')
        l += 1
        add.write('cross ; ' + ' '.join(map(str, [ssnear[n], 0])) + ' ; ' +
                  ' '.join(map(str, [a_sc - 1, 0])) + '\n')

# tn.close()
outpdb[0][3] += 'T'

for i, c in enumerate(atomsh):
    if i:
        allPath[c[0]][c[1]] = i
for vh in allPath:
    k = revers(vh, 0)
    if k == -1:
        allPath[vh] = allPath[vh][::-1]
    for i, c in enumerate(allPath[vh]):
        if c not in ['-', 'd', 's', 'DELE']:
            ntmp = c
            ctmp = 0
            allPath[vh][i] = (c, 0)
        elif c in ['d', 's']:
            ctmp += 1
            allPath[vh][i] = (ntmp, ctmp)
    if k == -1:
        allPath[vh] = allPath[vh][::-1]


for i in add_ends:
    p_i = allPath[i[1]][i[2]]
    if p_i == 'DELE':
        p_i = allPath[i[1]][i[2] - revers(i[1], 0)]
    add.write(i[0] + ' ; ' + str(p_i[0]) + ' ' + str(p_i[1]) + '\n')


def nsearch(b, a, direct):
    # determine direct
    i = allPath[b][a]
    while i == 'DELE':
        i = nsearch(b, a + direct, direct)
    return i


def wr_restr(i1, i2):
    global l, R_DONE
    if i1[1] not in [0, 1, HC - 1] or i2[1] not in [0, 1, HC - 1]:
        raise Exception('ADMIN: PY1. Restraints error')
    if (i1[0], i2[0]) not in R_DONE and (i2[0], i1[0]) not in R_DONE:
        # print i1, i2, l, len(R_DONE)
        cnct.write(str(i1[0]) + '\t' + str(i2[0]) +
                   '\t1\t' + str(l) + '\t1\t1.8\t1.85\t1.9\t1.4\n')
        l += 1
        R_DONE.append((i1[0], i2[0]))
        prconnect(i1[0], i2[0])


def wr_add(i1, i2):
    global T_DONE
    if (i1, i2) not in T_DONE and (i1, i2) not in T_DONE:
        add.write(
            'cross ; ' + ' '.join(map(str, i1)) + ' ' +
            ' '.join(map(str, i2)) + '\n')
        T_DONE.append((i1, i2))

PAIR_DONE = []
R_DONE = []
T_DONE = []

for pair in trans:
    [[b1, b2], [a1, a2]] = pair
    if not [[b2, b1], [a2, a1]] in PAIR_DONE:
        PAIR_DONE.append(pair)
        if not (a1 % HC or a2 % HC):
            i1 = nsearch(b1, a1, 1)
            i2 = nsearch(b2, a2, 1)
            wr_add(i1, i2)
            wr_restr(i1, i2)
            if [[b1, b2], [a1 - 1, a2 - 1]] in trans:
                i1 = nsearch(b1, a1 - 1, -1)
                i2 = nsearch(b2, a2 - 1, -1)
                wr_add(i1, i2)
        elif a1 % HC == HC - 1 and a2 % HC == HC - 1:
            if [[b1, b2], [a1 + 1, a2 + 1]] not in trans:
                i1 = nsearch(b1, a1, -1)
                i2 = nsearch(b2, a2, -1)
                wr_add(i1, i2)
#                print b1, a1, '|', b2, a2
#                print allPath[b1][a1-2:a1+3]
#                print i1, i2
                wr_restr(i1, i2)
        else:
            i1 = nsearch(b1, a1, 1)
            i2 = nsearch(b2, a2, 1)
            wr_add(i1, i2)
            wr_restr(i1, i2)

# print R_DONE, len(R_DONE)

# if a1 % HC == HC - 1 and a2 % HC == HC - 1: #and not pair
#        direct = -1
#        i1 = nsearch(b1, a1, direct)
#        i2 = nsearch(b2, a2, direct)
#        print 'add', i1, i2


#    i1 = [i for i, x in enumerate(atomsh) if x == (b1, a1)]
#    i2 = [i for i, x in enumerate(atomsh) if x == (b2, a2)]
#    print b1, a1, i1, '|', b2, a2, i2
#    if not (a1 % HC or a2 % HC):
#        print a1, a2
#        print i1, i2
#        i1 = nsearch(b1, a1, 1, 0)
#        i2 = nsearch(b2, a2, 1, 0)
#        print 'add ', i1, 0, i2, 0
#        i3= nsearch(b1, a1 - 1, -1, -1)
#        i4= nsearch(b2, a2 - 1, -1, -1)
#        print 'add ', i3, 0, i4, 0


#        print 'add ',

#    elif a1 % HC == HC - 1 and a2 % HC == HC - 1:

pdb = open(args.output, 'w')
prevB = False
for atom in outpdb:
    if len(atom) > 1:
        if (atom[3][0] == 'T' and
            atom[1] < a_sc - 1 and
            atom[1] > 1 and
            len(outpdb[atom[1]]) > 1 and
                not prevB):

            if outpdb[atom[1]][3] != 'S' and outpdb[atom[1] - 2][3] != 'S':
                new_name = TtoB(atom[3], atom[1])
                prevB = True
                outpdb[atom[1] - 1][3] = new_name
                outmap[atom[1] - 1][1] = new_name
                if new_name[0] == 'B':
                    if outpdb[atom[1]][3] == 'T':
                        outpdb[atom[1]][3] = 'B'
                        outmap[atom[1]][1] = 'B'
                    elif outpdb[atom[1]][3][0] == 'T':
                        outpdb[atom[1]][3] = 'B' + outpdb[atom[1]][3][1:]
                        outmap[atom[1]][1] = 'B' + outmap[atom[1]][1][1:]
        elif prevB:
            prevB = False
    outatom(atom)
for atom in outmap:
    atom[3] = vhNums.index(atom[3])
    mapf.write("{0[0]:s} {0[1]:s}:{0[2]:d}:{0[3]:d}:{0[4]:d}\n".format(atom))

# for i in range(len(connst)):
#    add.write('cross ; ' + ' '.join(map(str, connst[i])) + ' ; ' +
#              ' '.join(map(str, needst[i])) + '\n')
cnct.write('; scaffold crossovers\n')
if not endtmp:
    cnct.write('1\t' + str(num - 1) + '\t1\t' + str(l) +
               '\t1\t0.33\t0.35\t0.37\t2.5\t; 5\' - 3\' ends\n')
l += 1

for i in range(len(scconn)):
    if not (isinstance(scconn[i], int) and isinstance(scneed[i], int)):
        raise Exception('ADMIN: PY1. Error with scaffold restraints')
    cnct.write(
        str(scconn[i]) + '\t' + str(scneed[i]) + '\t1\t' + str(l) +
        '\t1\t0.33\t0.35\t0.37\t2.5\n'
    )
    l += 1
    add.write('scaf ; ' + str(scconn[i]) + ' ; ' + str(scneed[i]) + '\n')
cnct.write('\n')

add.close()

cnct.write('; lattice restraints\n')
if len(connects) != len(needed):
    print('Smth with connects is wrong')  # ~ error
cnct.write('\n')

template = "{0[0]:<6s}{0[1]:>5d}{0[2]:>5d}"
for i in range(1, num):     # doesn't work for ss oligs
    outpdbc.append(template.format(['CONECT', i, i + 1]) + '\n')

for conect in outpdbc:
    pdb.write(conect)

pdb.close()
mapf.close()

# ---------- lattice restraints ---------

name = 'H'

if SQ:
    [a, b, c] = ['4.34', '4.4', '4.43']
else:
    [a, b, c] = ['3.64', '3.81', '4.16']

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
                cnct.write(
                    str(triples[x][i]) + '\t' + str(triples[x][pair[1]]) +
                    '\t1\t' + str(l) + '\t1\t' + a + '\t' + b + '\t' + c +
                    '\t1.5\n')
                l += 1
cnct.close()
