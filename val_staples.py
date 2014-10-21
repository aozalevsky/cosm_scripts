import sys
import re

def val_staples(csv):
    for line in csv:
        line = line.split(',')
        print line
        for i in [0,1]:
            if not point.match(line[i]):
                return False
        if not seq.match(line[2]):
            return False
        if not leng.match(line[3]):
            return False
        if not color.match(line[4]):
            return False
    return True

point = re.compile("^[0-9]+\[[0-9]+\]$")
seq = re.compile("^[ATGC?]+$")
leng = re.compile("^[0-9]+$")
color = re.compile("^#.+$")

with open(sys.argv[1], 'r') as csv:
    print val_staples(csv)