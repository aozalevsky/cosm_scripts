#!/usr/bin/env python
# encoding: utf-8

import json
import sys


def load_json(fname):
    """ Load cadnano json file"""
    try:
        with open(fname, 'r') as f:
            obj = json.load(f)

    except IOError as e:
        raise Exception("ADMIN: %s" % e)

    return obj


def check_json_version(obj):
    """Check format of json. Currently we support only <2.0"""

    support = True

    if "format" in obj:
        if float(obj['format']) > 2.0:
            support = False

    if support is False:
        raise Exception('USER: Unsupported cadnano json version')


def val_lattice(f):
    numBases = len(obj['vstrands'][0]['scaf'])
    if numBases % 32 == 0 and numBases % 21 != 0:
        return 's'
    elif numBases % 21 == 0 and numBases % 32 != 0:
        return 'h'
    else:
        return 'u'  # unknown


if __name__ == '__main__':

    obj = load_json(sys.argv[1])

    try:
        print(val_lattice(obj))
    except:
        raise Exception('USER: Error in input json file')
