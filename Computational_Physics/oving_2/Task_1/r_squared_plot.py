import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

#data_100 = np.loadtxt('res_hon/p100.csv', delimiter='\n')[40:-40]
data_200 = np.loadtxt('res_hon/p200.csv', delimiter='\n')[40:-40]
data_300 = np.loadtxt('res_hon/p300.csv', delimiter='\n')[40:-40]
data_400 = np.loadtxt('res_hon/p400.csv', delimiter='\n')[40:-40]
data_600 = np.loadtxt('res_hon/p600.csv', delimiter='\n')[40:-40]
data_500 = np.loadtxt('res_hon/p500.csv', delimiter='\n')[40:-40]
data_700 = np.loadtxt('res_hon/p700.csv', delimiter='\n')[40:-40]
data_800 = np.loadtxt('res_hon/p800.csv', delimiter='\n')[40:-40]
data_900 = np.loadtxt('res_hon/p900.csv', delimiter='\n')[40:-40]
data_1000= np.loadtxt('res_hon/p1000.csv', delimiter='\n')[40:-40]
data_1500= np.loadtxt('res_hon/p1500.csv', delimiter='\n')[40:-40]

def get_r_value(p_inf_value):
    p_inf_values = np.array([data_200[p_inf_value], data_300[p_inf_value], data_400[p_inf_value], data_500[p_inf_value], data_600[p_inf_value], data_700[p_inf_value], data_800[p_inf_value], data_900[p_inf_value], data_1000[p_inf_value], data_1500[p_inf_value]])
    n_sites = np.array([200,300,400,500,600,700, 800, 900, 1000, 1500])

    slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(p_inf_values), np.log(n_sites))
    return rvalue**2

def calculate_r_squared():
    p = np.linspace(0.0,1.0,len(data_1000))
    r_squared_values = np.zeros(len(p))

    for i in range(0,len(p)):
        r_squared_values[i] = get_r_value(i)
    return p, r_squared_values


p, r_squared_values = calculate_r_squared()

start = 6450
end = 9000

#Get pc
pc = start + np.argmax(r_squared_values[start:end])
print "Probability that makes R*R larger: ", float(pc)/(len(data_1000))

plt.figure(1)
plt.plot(p[start:end], r_squared_values[start:end])
plt.show()
