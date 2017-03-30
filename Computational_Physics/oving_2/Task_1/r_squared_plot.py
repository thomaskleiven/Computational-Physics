import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

data_100 = np.loadtxt('results/p100.csv', delimiter='\n')[10:-10]
data_200 = np.loadtxt('results/p200.csv', delimiter='\n')[10:-10]
data_300 = np.loadtxt('results/p300.csv', delimiter='\n')[10:-10]
data_400 = np.loadtxt('results/p400.csv', delimiter='\n')[10:-10]
data_600 = np.loadtxt('results/p600.csv', delimiter='\n')[10:-10]
data_500 = np.loadtxt('results/p500.csv', delimiter='\n')[10:-10]
data_700 = np.loadtxt('results/p700.csv', delimiter='\n')[10:-10]
data_800 = np.loadtxt('results/p800.csv', delimiter='\n')[10:-10]
data_900 = np.loadtxt('results/p900.csv', delimiter='\n')[10:-10]
data_1000= np.loadtxt('results/p1000.csv', delimiter='\n')[10:-10]

def get_r_value(p_inf_value):
    p_inf_values = np.array([data_100[p_inf_value], data_200[p_inf_value], data_300[p_inf_value], data_400[p_inf_value], data_500[p_inf_value], data_600[p_inf_value], data_700[p_inf_value], data_800[p_inf_value], data_900[p_inf_value], data_1000[p_inf_value]])
    n_sites = np.sqrt(np.array([100**2, 200**2, 300**2, 400**2, 500**2, 600**2, 700**2, 800**2, 900**2, 1000**2]))

    slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(p_inf_values), np.log(n_sites))
    return rvalue**2

def calculate_r_squared():
    p = np.linspace(0.0,1.0,len(data_100))
    r_squared_values = np.zeros(len(p))

    for i in range(0,len(p)):
        r_squared_values[i] = get_r_value(i)
    return p, r_squared_values


p, r_squared_values = calculate_r_squared()

start = 4900
end = 9000

#Get pc
pc = start + np.argmax(r_squared_values[start:end])
print "Probability that makes R*R larger: ", float(pc)/(len(data_100))

plt.figure(1)
plt.plot(p[start:end], r_squared_values[start:end])
plt.show()
