import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

data_300 = np.loadtxt('average_sq_p/p300__averaged.csv', delimiter='\n')[10:-10]
data_500 = np.loadtxt('average_sq_p/p500__averaged.csv', delimiter='\n')[10:-10]
data_700 = np.loadtxt('average_sq_p/p700__averaged.csv', delimiter='\n')[10:-10]
data_800 = np.loadtxt('average_sq_p/p800__averaged.csv', delimiter='\n')[10:-10]
data_1000 = np.loadtxt('average_sq_p/p1000__averaged.csv', delimiter='\n')[10:-10]
data_2000 = np.loadtxt('average_sq_p/p2000__averaged.csv', delimiter='\n')[10:-10]


start = 1
p_c = 0.499932975871

n_sites = np.array([300,500,700,800, 1000,2000])
p_inf_values = np.array([data_300[p_c*10000], data_500[p_c*10000], data_700[p_c*10000], data_800[p_c*10000], data_1000[p_c*10000], data_2000[p_c*10000]])


#Fit line to p_inf_values
slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(n_sites[start:]), np.log(p_inf_values[start:]))

print "Slope of p_inf as a function of correlation length: ", slope

plt.figure(2)
plt.plot(n_sites, np.exp(intercept)*n_sites**slope)
plt.plot(n_sites, p_inf_values, "x")
plt.xscale("log")
plt.yscale("log")
plt.show()
