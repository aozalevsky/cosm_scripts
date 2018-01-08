#!/usr/bin/env python

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import sys

matplotlib.use('pdf')


f = open(sys.argv[1])
lines = f.readlines()
f.close()

text = []
verts = []
dele = []
ins = []
codes = [Path.MOVETO]
lenx = 0
leny = 0
for l in lines:
    l = l.split()
    if l[0] == 'L':
        x, y = l[1].split(":")
        lenx = max(float(x), lenx)
        leny = max(float(y), leny)
        verts.append((float(x), float(y)))
        codes.append(Path.LINETO)
    elif l[0] == 'T':
        t, t1, y, x = l[1].split(":")
        text.append([float(x), float(y), t, t1])
    elif l[0] == 'D':
        x, y = l[1].split(":")
        dele.append([float(x), float(y)])
    elif l[0] == 'I':
        t, t1, y, x = l[1].split(":")
        ins.append([float(x), float(y), t, t1])
    else:
        raise Exception("Error map file")

del codes[-1]
path = Path(verts, codes)

fig = plt.figure(figsize=(lenx / 5.0, leny / 2.0))
ax = plt.subplot(111)

# ax.annotate("",xy=(0.2, 0.2), xycoords='data',xytext=(0.8, 0.8),
# textcoords='data',arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),)

ax = fig.add_subplot(111)
patch = patches.PathPatch(
    path, facecolor='none', edgecolor='orange', alpha=0.5, lw=1)

ax.add_patch(patch)

bbox_props = dict(boxstyle="round,pad=0.3",
                  ec="gray", lw=0.5, facecolor="gray", alpha=0.15)
bbox_small = dict(boxstyle="round,pad=0.3",
                  ec="gray", lw=0.5, facecolor="gray", alpha=0.15)

for la in text:
    ax.text(la[0], la[1], la[2], fontsize=6, bbox=bbox_props, ha='center')
    ax.text(la[0], la[1] - 0.20, la[3], fontsize=3, ha='center')
for la in dele:
    x = la[0]
    y = la[1]
    ax.annotate("",  xy=(x - 0.4, y - 0.15), xycoords='data', xytext=(
        x + 0.4, y + 0.15), textcoords='data', arrowprops=dict(arrowstyle="-"))
    ax.annotate("",  xy=(x + 0.4, y - 0.15), xycoords='data', xytext=(
        x - 0.4, y + 0.15), textcoords='data', arrowprops=dict(arrowstyle="-"))
for la in ins:
    x = la[0]
    y = la[1]
    ax.annotate("",  xy=(x - 0.4, y + 0.4), xycoords='data', xytext=(
        x, y + 0.2), textcoords='data', arrowprops=dict(arrowstyle="-"))
    ax.annotate("",  xy=(x + 0.4, y + 0.4), xycoords='data', xytext=(
        x, y + 0.2), textcoords='data', arrowprops=dict(arrowstyle="-"))
    ax.text(la[0] - 0.2, la[1] + 0.4, la[2],
            fontsize=6, bbox=bbox_props, ha='center')
    ax.text(la[0] - 0.2, la[1] + 0.6, la[3], fontsize=3, ha='center')

ax.relim()
ax.autoscale_view(None, True, True)
plt.axis('off')
plt.savefig(sys.argv[2])
# ,bbox_inches='tight')
