import numpy as np
import matplotlib as mpl
import sys
from os import listdir
from os.path import isfile, join
mpl.rcParams['svg.fonttype'] = "none"
mpl.rcParams['font.size'] = 18
from matplotlib import pyplot as plt

#colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"]
colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"]
path = sys.argv[1]

files = [f for f in listdir(path) if isfile(join(path, f))]


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
count = 0
for filename in files:
    data = np.loadtxt(path+"/"+filename)
    x = np.linspace(0.0, 1.0, len(data))
    ax.plot(x, data, color=colors[count%len(colors)], label=filename)
    count += 1


ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")
ax.set_xlabel("\$p\$")
ax.set_ylabel("\$P_\infty\$")
ax.legend(loc="lower right", frameon=False, labelspacing=0.05)
plt.show()
