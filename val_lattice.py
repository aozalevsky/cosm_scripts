#!/usr/bin/env python
# encoding: utf-8

import json
import sys


def val_lattice(f):
    obj = json.load(f)
    numBases = len(obj['vstrands'][0]['scaf'])
    if numBases % 32 == 0 and numBases % 21 != 0:
        return 's'
    elif numBases % 21 == 0 and numBases % 32 != 0:
        return 'h'
    else:
        return 'u'  # unknown

if __name__ == '__main__':

    try:
        with open(sys.argv[1], 'r') as f:
            print(val_lattice(f))
    except:
        raise Exception('USER: Error in input json file')
