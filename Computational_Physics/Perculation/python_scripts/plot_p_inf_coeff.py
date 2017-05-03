import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

import matplotlib as mpl
mpl.rcParams['svg.fonttype'] = "none"
mpl.rcParams['font.size'] = 18
mpl.rcParams['legend.numpoints'] = 1

data_300 = np.loadtxt('square/p300__averaged.csv', delimiter='\n')[40:-40]
data_500 = np.loadtxt('square/p500__averaged.csv', delimiter='\n')[40:-40]
data_700 = np.loadtxt('square/p700__averaged.csv', delimiter='\n')[40:-40]
data_800 = np.loadtxt('square/p800__averaged.csv', delimiter='\n')[40:-40]
data_1000 = np.loadtxt('square/p1000__averaged.csv', delimiter='\n')[40:-40]
data_2000 = np.loadtxt('square/p2000__averaged.csv', delimiter='\n')[40:-40]


start = 1
#p_c = 0.499932975871
#p_c = 0.653619302949
#p_c = 0.346481233244

#tri = 10339
#square = 14919
#hon = 19503


n_sites = np.array([300,500,700,800, 1000,2000])
p_inf_values = np.array([data_300[14919], data_500[14919], data_700[14919], data_800[14919], data_1000[14919], data_2000[14919]])


#Fit line to p_inf_values
slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(n_sites[start:]), np.log(p_inf_values[start:]))

print ( "Slope of p_inf as a function of correlation length: ", slope)

a, V = np.polyfit(np.log(n_sites[start:]), np.log(p_inf_values[start:]), 1, cov=True)

print (a)
print (V)


fig = plt.figure(2)
ax = fig.add_subplot(1,1,1)
ax.plot(n_sites, np.exp(intercept)*n_sites**slope, label="\$Fit a \cdot x+b$")
ax.plot(n_sites, p_inf_values, "x", label="\$data$")

ax.set_yscale('log')
ax.set_xscale('log')
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")
ax.set_xlabel("\$log xi $")
ax.set_ylabel("\$log p_{\infty}$")
ax.legend(loc="upper right", frameon=False, labelspacing=0.05)

plt.show()
