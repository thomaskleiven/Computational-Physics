from matplotlib import pyplot as plt
import numpy as np

import matplotlib as mpl
mpl.rcParams['svg.fonttype'] = "none"
mpl.rcParams['font.size'] = 18

mode = np.loadtxt("../crankNicolsonScheme.csv", delimiter=",")
mode = np.transpose(mode)
x = np.linspace(0.0,1.0, len(mode[0]))

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")
ax.set_xlabel("\$ x (a.u) $")
ax.set_ylabel("\$\psi (a.u)$")

ax.plot(x, mode[12], color="black", label="\$Numerical$")
ax.legend(loc="upper right", frameon=False, labelspacing=0.05)

plt.show()
