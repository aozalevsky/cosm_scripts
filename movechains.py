#!/usr/bin/python

import sys

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')

call = list()

for line in infile:
    if line[:4] == 'ATOM':
        beg = line[:29]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        end = line[55:]
        if (x, y, z) not in call:
            call.append((x, y, z))
            outfile.write(line)
        else:
            while (x, y, z) in call:
                x += 0.2
            outfile.write(
                " ".join(
                    [
                        beg,
                        "{0[0]:>8.3f}{0[1]:>8.3f}{0[2]:>8.3f}".format(
                            (x, y, z)),
                        end
                    ]))
            call.append((x, y, z))
    elif line.strip():
        outfile.write(line)

infile.close()
outfile.close()
