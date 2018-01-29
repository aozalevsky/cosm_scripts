#!/usr/bin/env python
# encoding: utf-8

import sys
import json
import re
import pyparsing as pg
from collections import Counter


def load_json(fname):
    """ Load cadnano json file"""
    try:
        with open(fname, 'r') as f:
            obj = json.load(f)

    except IOError as e:
        raise Exception("ADMIN: %s" % e)

    check_json_version(obj)

    return obj


def check_json_version(obj):
    """Check format of json. Currently we support only <2.0"""

    support = True

    if "format" in obj:
        if float(obj['format']) > 2.0:
            support = False

    # Check that it's really v2 format
    if 'vstrads' in obj:
        pass

    if support is False:
        raise Exception('USER: Unsupported cadnano json version')


def get_lattice_type(obj):
    numBases = len(obj['vstrands'][0]['scaf'])
    if numBases % 32 == 0 and numBases % 21 != 0:
        return 's'
    elif numBases % 21 == 0 and numBases % 32 != 0:
        return 'h'
    else:
        raise Exception(
            'USER: lattice neither honeycomb nor square')
        # return 'u'  # unknown


def find_duplicates(arr):
    return [item for item, count in Counter(arr).items() if count > 1]


lattice_types = [
    'h', 'hex', 'hexagonal',
    's', 'sq', 'square',
]


class COSMError(Exception):

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


compl = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
    'a': 'T',
    't': 'A',
    'g': 'C',
    'c': 'G'
}


def parse_staples(fn):
    try:
        with open(fn, 'r') as f:
            raw = f.readlines()
    except:
        raise COSMError('USER: Error in sequence input file')

    # example coordinate: 1[23]
    coord = pg.Word(pg.nums).setName('coordinates in format XX[YY].')
    coord.setParseAction(lambda s, l, t: [int(t[0])])  # Convert num to int
    # example sequence: AGGTTAGATCG
    seq = pg.Regex("[ATGC]+").setName('DNA sequence composed of ATGC.')
    # example length: 8
    length = pg.Word(pg.nums).setName('integer staple length.')
    # example color: #54aa8d
    color = pg.Regex('#[0-9A-z]{6}').setName('color in format #AABBCC.')

    ob = pg.Suppress('[')
    cb = pg.Suppress(']')
    com = pg.Suppress(',')

    # example line
    # 0[116],1[107],TGTCTCAGCTGCATCGCAAGACATCATCAAAGG,33,#170fde
    validLine = (
        coord + ob + coord + cb + com +
        coord + ob + coord + cb + com +
        seq + com +
        length + com +
        color +
        pg.LineEnd()
    )

    oligs = list()

    for i in range(len(raw)):
        line = raw[i]
        try:
            oligs.append(validLine.parseString(line)[:-1])  # skip line end
        except pg.ParseException as pe:
            # Skip header line
            if i == 0 and re.match('Start', line):
                continue
            else:
                msg, error = pe.markInputline(
                    '?'), 'Invalid staples. ' + pe.msg
                raise COSMError(
                    'USER: %s. %s' % (msg, error))

    return oligs


def parse_sequence(seqfn):
    try:
        with open(seqfn, 'r') as f:
            raw = f.readlines()
    except:
        raise COSMError('USER: Error in sequence input file')

    reg_seq = pg.Regex("[ATGC]+").setName('DNA sequence composed of ATGC.')

    validLine = reg_seq + pg.LineEnd()

    seq = ''

    for line in raw:
        try:
            seq += validLine.parseString(line.upper())[0]
        except pg.ParseException as pe:
            error = 'Illegal symbols in DNA sequence.'
            msg = pe.markInputline('?')
            raise COSMError(
                'USER: %s. %s' % (msg, error))

    return seq


def pretty(t, flag, msg=None, error=None):
    print('-' * 80)
    print('TYPE:', t)
    print('MSG:', msg)
    print('ERROR:', error)
    print('-' * 80)


if __name__ == '__main__':

    obj = load_json(sys.argv[1])

    print(get_lattice_type(obj))
