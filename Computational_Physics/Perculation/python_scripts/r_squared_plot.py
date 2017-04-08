import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

data_300 = np.loadtxt('average_sq_p/p300__averaged.csv', delimiter='\n')[40:-40]
data_500 = np.loadtxt('average_sq_p/p500__averaged.csv', delimiter='\n')[40:-40]
data_700 = np.loadtxt('average_sq_p/p700__averaged.csv', delimiter='\n')[40:-40]
data_800 = np.loadtxt('average_sq_p/p800__averaged.csv', delimiter='\n')[40:-40]
data_1000= np.loadtxt('average_sq_p/p1000__averaged.csv', delimiter='\n')[40:-40]
data_2000= np.loadtxt('average_sq_p/p2000__averaged.csv', delimiter='\n')[40:-40]

def get_r_value(p_inf_value):
    p_inf_values = np.array([data_300[p_inf_value], data_500[p_inf_value], data_700[p_inf_value], data_800[p_inf_value], data_1000[p_inf_value], data_2000[p_inf_value]])
    n_sites = np.array([300, 500 ,700, 800, 1000, 2000])

    slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(p_inf_values), np.log(n_sites))
    return rvalue**2

def calculate_r_squared():
    p = np.linspace(0.0,1.0,len(data_1000))
    r_squared_values = np.zeros(len(p))

    for i in range(0,len(p)):
        r_squared_values[i] = get_r_value(i)
    return p, r_squared_values


p, r_squared_values = calculate_r_squared()

start = 14800
end = 24000

#Get pc
pc = start + np.argmax(r_squared_values[start:])
print "Probability that makes R*R larger: ", float(pc)/(len(data_1000))

plt.figure(1)
plt.plot(p[start:end], r_squared_values[start:end])
plt.show()
