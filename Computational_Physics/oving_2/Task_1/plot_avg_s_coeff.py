import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

start = 0
end = 10
p_c = 0.499932975871
n_sites = np.array([300,500,700,800,1000,2000])

data_300 = np.loadtxt('average_sq_avg/avg300__averaged.csv', delimiter='\n')[40:-40]
data_500 = np.loadtxt('average_sq_avg/avg500__averaged.csv', delimiter='\n')[40:-40]
data_700 = np.loadtxt('average_sq_avg/avg700__averaged.csv', delimiter='\n')[40:-40]
data_800 = np.loadtxt('average_sq_avg/avg800__averaged.csv', delimiter='\n')[40:-40]
data_1000 = np.loadtxt('average_sq_avg/avg1000__averaged.csv', delimiter='\n')[40:-40]
data_2000 = np.loadtxt('average_sq_avg/avg2000__averaged.csv', delimiter='\n')[40:-40]

#Calculate p_max
def calculate_p_max():
    values = np.array([data_300, data_500, data_700, data_800, data_1000, data_2000])
    max_s_vals = np.zeros(len(values))
    p_max = np.zeros(len(values))
    for i in range(0, len(values)):
        max_s_vals[i] = np.max(values[i])
        p_max[i] = np.argmax(values[i])

    return p_max/(len(data_800)), max_s_vals

pmax, max_s_vals = calculate_p_max()

#Plot abs(p_max - p_c) against correlation lenght
def plot_ny():
    p_max, max_s_vals = calculate_p_max()
    p_max_p_c = abs(p_max - p_c)

    #Fit line
    slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(n_sites[start:end]), np.log(p_max_p_c[start:end]))
    print "Slope of |pmax-pc| as a function of correlation length: ", slope

    #Fit line
    a, V = np.polyfit(np.log(n_sites[start:end]), np.log(p_max_p_c[start:end]), 1, cov=True)

    print a
    print V

    plt.figure(10)
    plt.plot(n_sites, np.exp(intercept)*n_sites**slope)
    plt.plot(n_sites, p_max_p_c, 'x')
    plt.xscale("log")
    plt.yscale("log")
    plt.show()

#Fit line to max(<s>)
slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(n_sites[start:end]), np.log(max_s_vals[start:end]))
print "Slope of max(<s>) as a function of correlation length: ", slope

#a, V = np.polyfit(np.log(n_sites[start:end]), np.log(max_s_vals[start:end]), 1, cov=True)

plot_ny()
