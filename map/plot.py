import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import sys


f = open(sys.argv[1])
lines = f.readlines()
f.close()

verts=[]
codes=[Path.MOVETO]
for l in lines:
    l.strip()
    x,y=l.split(":")
    verts.append((float(x),float(y)))
    codes.append(Path.LINETO)

f = open(sys.argv[2])
lines = f.readlines()
f.close()
text=[]
for l in lines:
    l.strip()
    t,t1,x,y=l.split(":")
    text.append([float(x),float(y),t])

del codes[-1]
print len(verts)
print len(codes)
path = Path(verts, codes)

fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(111)
patch = patches.PathPatch(path, facecolor='none',edgecolor='orange', alpha=0.5, lw=1)
#ax.set_xlim('auto')
#ax.set_ylim('auto')

ax.add_patch(patch)

for la in text:
    ax.text(la[0], la[1], la[2], fontsize=6, bbox={'facecolor':'gray', 'alpha':0.1,
        'pad':2})

ax.relim()
ax.autoscale_view(None,True,True)
plt.axis('off')
plt.savefig(sys.argv[3],bbox_inches='tight')
