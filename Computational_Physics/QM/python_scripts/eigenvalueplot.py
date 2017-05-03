from matplotlib import pyplot as plt
import numpy as np

import matplotlib as mpl
mpl.rcParams['svg.fonttype'] = "none"
mpl.rcParams['font.size'] = 18

data = np.loadtxt("../eigenvalues/eigenvalues_100.csv", delimiter=",")
x = np.zeros(len(data))

for n in range(0,len(data)):
    x[n] = (n*np.pi)**2

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")
ax.set_xlabel("\$log xi $")
ax.set_ylabel("\$log p_{\infty}$")

ax.plot(x, x, "+", color="red", label="\$Exact$")
ax.plot(x, data, "x", color="black", label="\$Numerical$")
ax.legend(loc="upper left", frameon=False, labelspacing=0.05)

plt.show()
