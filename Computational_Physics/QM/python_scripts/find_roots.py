from matplotlib import pyplot as plt
import numpy as np
from scipy import optimize, fmin
import matplotlib as mpl
mpl.rcParams['svg.fonttype'] = "none"
mpl.rcParams['font.size'] = 18

v0 = 1000

#Callable function
def eigvecs(x):
    k = np.sqrt(x)
    kappa = np.sqrt(v0-x)
    first = np.exp(kappa/3)
    second = (kappa*np.sin(k/3) + k*np.cos(k/3))**2
    third = np.exp(-kappa/3)
    fourth = (kappa*np.sin(k/3) - k*np.cos(k/3))**2

    return first*second - third*fourth

#Find minimum
xopt = optimize.fmin(eigvecs, 600) #Function, initial guess
print "Optimal x-value: ", xopt[0]

#Find root
root = optimize.brentq(eigvecs, 648.507041931, 600) #Function, inital guess
print "Root: ", root


a = np.zeros(999)
for i in range(0,999):
    a[i] = eigvecs(i)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.axhline(0, color="black")
ax.plot(a)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")
ax.set_xlabel("\$ x (a.u) $")
ax.set_ylabel("\$ \lambda (a.u)$")
#ax.legend(loc="upper right", frameon=False, labelspacing=0.05)
plt.show()
