#!/usr/bin/env python
# encoding: utf-8

import json, sys


def val_lattice(file):
    lines = file.readlines()
    stringl = ""
    for line in lines:
        stringl += line
    obj = json.loads(stringl)
    numBases = len(obj['vstrands'][0]['scaf'])
    if numBases % 32 == 0:
        return 's'
    elif numBases % 21 == 0:
        return 'h'
    else:
        return None

with open(sys.argv[1], 'r') as file:
    print val_lattice(file)