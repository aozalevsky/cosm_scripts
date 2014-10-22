#!/usr/bin/env python
# encoding: utf-8


def val_lattice(f):
    import json
    obj = json.load(f)
    numBases = len(obj['vstrands'][0]['scaf'])
    if numBases % 32 == 0:
        return 's'
    elif numBases % 21 == 0:
        return 'h'
    else:
        return None
if __name__ == '__main__':
    import sys
    with open(sys.argv[1], 'r') as f:
        print val_lattice(f)
