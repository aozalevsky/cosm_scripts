#!/usr/bin/env python
# encoding: utf-8

import re
import sys

def s6(line, old):
    n = None
    c1 = 100500
    i = 0
    old2 = []
    while (not n or c1 in old) and i <= 3:
        n = re.search(r'[T, B]6A', line)
        if n:
            c1 = n.start()
            line = line[:c1 + 1] + '7' + line[c1 + 2:]
            if c1 in old:
                old2.append(c1)
#        print c1, old, n, i
        i += 1
    for i in old2:
        line = line[:i + 1] + '6' + line[i + 2:]
    if n:
        print line.strip()
        return line, c1
    return [None, None]

ccg7 = open(sys.argv[1], 'r')
for line in ccg7:
    i = 0
    k = True
    while i <= 3 and k:
        n = re.search(r'[T, B]7A', line)
        if n:
            c1 = n.start()
            line = line[:c1 + 1] + '8' + line[c1 + 2:]
            i += 1
            k = True
        else: k = False
    print line.strip()
    [nline, cn] = s6(line, [])  # 1
    old = [cn]
    if nline:
        [n2line, cn] = s6(nline, old)  # 2
        old.append(cn)
        if n2line:
            [n3line, cn] = s6(n2line, old)  # 3
            old.append(cn)
    if len(old) > 1:
        if old[1] and not old[2]:
#            print line, old[0]
            s6(line, [old[0]])
    if len(old) > 2:
        if old[2]:
            s6(line, old[0:2])
            s6(line, [old[0]])
            s6(nline, [old[1]])
            n = re.search(r'[B, T]7A', n3line)
            c1 = n.start()
            line = n3line[:c1 + 1] + '6' + n3line[c1 + 2:]
            print line.strip()
#    else:
#        print line

ccg7.close()
