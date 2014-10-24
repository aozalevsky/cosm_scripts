#!/usr/bin/env python
# encoding: utf-8


def validate_json(f):
    # initiate variables
    flag, msg, error = True, 'JSON is OK.', ''

# ------------- .json analysis ----------------
    import json

    try:
        obj = json.load(f)
        strands = obj["vstrands"]

# create dictionaries (keyed by virtual helix #) of
# row/col, scaf array, stap array

        vhToScaf = {}
        vhToStap = {}
        vhNums = []
        Rows = {}
        Cols = {}
        trans = []
        trans_sc = []
        for strand in strands:
            num = strand["num"]
            vhNums.append(num)
            Rows[num] = strand["row"]
            Cols[num] = strand["col"]
            scaf = strand["scaf"]
            stap = strand["stap"]
            vhToScaf[num] = scaf
            vhToStap[num] = stap

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
                if (base[2] == -1) & (base[0] == vh):
                    scafPath.append('e3')

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
                raise RuntimeError
        for base in Path['scaf'][vh]:
            if base not in (range(0, 10) + ['-', 'e5', 'e3']):
                raise RuntimeError
        for base in Path['stap'][vh]:
            if base not in (range(0, 10) + ['-', 'e5', 'e3']):
                raise RuntimeError

    except:
        flag = False
        msg = ''
        error = 'JSON is invalid.'

    return (flag, msg, error)


def validate_lattice(f):
    # initiate variables
    flag, msg, error = True, '', ''

    flag, msg, error = validate_json(f)

    if flag:
        f.seek(0)

        import json
        obj = json.load(f)
        numBases = len(obj['vstrands'][0]['scaf'])

        if numBases % 32 == 0 and numBases % 21 != 0:
            msg = 's'
        elif numBases % 21 == 0 and numBases % 32 != 0:
            msg = 'h'
        else:
            flag = False
            error = 'Unable to detect lattice type.'

    return (flag, msg, error)


def validate_staples(f):
    # initiate variables
    flag, msg, error = True, 'Staples are OK', ''

    import pyparsing as pg
    # example coordinate: 1[23]
    coord = pg.Regex('\d+\[\d+\],').setName('coordinates in format XX[YY].')
    # example sequence: AGGTTAGATCG
    seq = pg.Regex("[ATGC]+,").setName('DNA sequence composed of ATGC.')
    # example length: 8
    length = pg.Regex('\d+,').setName('integer staple length.')
    # example color: #54aa8d
    color = pg.Regex('#[0-9A-z]{6}').setName('color in format #AABBCC.')

    # example line
    # 0[116],1[107],TGTCTCAGCTGCATCGCAAGACATCATCAAAGG,33,#170fde
    validLine = coord + coord + seq + length + color + pg.LineEnd()

    for line in f:
        try:
            validLine.parseString(line)
        except pg.ParseException, pe:
            flag, msg, error = False, pe.markInputline('?'), 'Invalid staples. ' + pe.msg
            break

    return (flag, msg, error)


def validate_sequence(f):
    flag, msg, error = True, 'Sequence is OK', ''

    import pyparsing as pg
    # example coordinate: 1[23]
    seq = pg.Regex("[ATGC]+").setName('DNA sequence composed of ATGC.')

    validLine = seq + pg.LineEnd()

    for line in f:
        try:
            validLine.parseString(line.upper())
        except pg.ParseException, pe:
            flag = False
            error = 'Illegal symbols in DNA sequence.'
            msg = pe.markInputline('?')

    return (flag, msg, error)


def pretty(t, flag, msg=None, error=None):
    print '-' * 80
    print 'TYPE:', t
    print 'MSG:', msg
    print 'ERROR:', error
    print '-' * 80


if __name__ == '__main__':
    import argparse as ag

    parser = ag.ArgumentParser()

    parser.add_argument('-j', '--json', type=ag.FileType('r'),
                        help='input cadnano json file')
    parser.add_argument('-l', '--lattice', type=ag.FileType('r'),
                        help='input cadnano json file')
    parser.add_argument('-s', '--sequence', type=ag.FileType('r'),
                        help='input cadnano json file')
    parser.add_argument('-st', '--staples', type=ag.FileType('r'),
                        help='input cadnano json file')
    args = parser.parse_args()

    if args.json:
#	print validate_json(args.json)
        pretty('JSON', *validate_json(args.json))
        args.json.close()
    if args.lattice:
        pretty('Lattice', *validate_lattice(args.lattice))
        args.lattice.close()
    if args.sequence:
        pretty('Sequence', *validate_sequence(args.sequence))
        args.sequence.close()
    if args.staples:
        pretty('Staples', *validate_staples(args.staples))
        args.staples.close()
