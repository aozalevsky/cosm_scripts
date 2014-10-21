#!/usr/bin/env python
# encoding: utf-8

#import __main__
#__main__.pymol_argv = ['pymol', '-qc']
#import pymol
import json
import sys

# ------------- .json analysis ----------------

def val_json(file):
    lines = file.readlines()
    stringl = ""
    try:
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
#        Path['scaf'][vh][5] = 'a'
        for base in allPath[vh]:
            if base not in ['-', 's', 'd', 'o']:
                return False
        for base in Path['scaf'][vh]:
            if base not in (range(0, 10) + ['-', 'e5', 'e3']):
                return False
        for base in Path['stap'][vh]:
            if base not in (range(0, 10) + ['-', 'e5', 'e3']):
                return False

        return True
    except:
        return False

with open(sys.argv[1], 'r') as file:
    print val_json(file)