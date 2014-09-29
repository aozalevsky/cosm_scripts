#!/usr/bin/env python
# encoding: utf-8

import sys

beg = False
end = False

template = "{0[0]:<6s}{0[1]:>5d}{0[2]:>5d}"

out = open(sys.argv[1] + '_end.pdb', 'a')

with open(sys.argv[1] + '/' + sys.argv[1] + '_r', 'r') as f:
    for line in f:
        line = line.split()
        if len(line) > 2:
            if line[1] == 'staple':
                beg = True
            elif line[1] == 'scaffold':
                end = True
            elif beg and not end:
                line[0] = int(line[0])
                line[1] = int(line[1])
                out.write(template.format(['CONECT'] + line[:2]) + '\n')
