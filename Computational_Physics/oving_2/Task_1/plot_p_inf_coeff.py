import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

data_100 = np.loadtxt('results/p100.csv', delimiter='\n')
data_200 = np.loadtxt('results/p200.csv', delimiter='\n')
data_300 = np.loadtxt('results/p300.csv', delimiter='\n')
data_400 = np.loadtxt('results/p400.csv', delimiter='\n')
data_500 = np.loadtxt('results/p500.csv', delimiter='\n')
data_600 = np.loadtxt('results/p600.csv', delimiter='\n')
data_700 = np.loadtxt('results/p700.csv', delimiter='\n')
data_800 = np.loadtxt('results/p800.csv', delimiter='\n')
data_900 = np.loadtxt('results/p900.csv', delimiter='\n')
data_1000 = np.loadtxt('results/p1000.csv', delimiter='\n')


start = 1

n_sites = np.sqrt(np.array([100**2, 200**2, 300**2, 400**2, 500**2, 600**2, 700**2, 800**2, 900**2, 1000**2]))
p_inf_values = np.array([data_100[4800], data_200[4800], data_300[4800], data_400[4800], data_500[4800], data_600[4800], data_700[4800], data_800[4800], data_900[4800], data_1000[4800]])


#Fit line to p_inf_values
slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(n_sites[start:]), np.log(p_inf_values[start:]))

print "Slope of p_inf as a function of correlation length: ", slope

plt.figure(2)
plt.plot(n_sites, np.exp(intercept)*n_sites**slope)
plt.plot(n_sites, p_inf_values, "x")
plt.xscale("log")
plt.yscale("log")
plt.show()
