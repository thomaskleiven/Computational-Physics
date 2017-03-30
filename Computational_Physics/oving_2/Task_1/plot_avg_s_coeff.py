import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

start = 0
end = 10
p_c = 0.499498997996
n_sites = np.sqrt(np.array([100**2, 200**2, 400**2, 500**2, 600**2, 700**2, 800**2, 1000**2]))

data_100 = np.loadtxt('results/avg100.csv', delimiter='\n')[10:-10]
data_200 = np.loadtxt('results/avg200.csv', delimiter='\n')[10:-10]
#data_300 = np.loadtxt('results/avg300.csv', delimiter='\n')[10:-10]
data_400 = np.loadtxt('results/avg400.csv', delimiter='\n')[10:-10]
data_500 = np.loadtxt('results/avg500.csv', delimiter='\n')[10:-10]
data_600 = np.loadtxt('results/avg600.csv', delimiter='\n')[10:-10]
data_700 = np.loadtxt('results/avg700.csv', delimiter='\n')[10:-10]
data_800 = np.loadtxt('results/avg800.csv', delimiter='\n')[10:-10]
#data_900 = np.loadtxt('results/avg900.csv', delimiter='\n')[10:-10]
data_1000= np.loadtxt('results/avg1000.csv', delimiter='\n')[10:-10]

#Calculate p_max
def calculate_p_max():
    values = np.array([data_100, data_200, data_400, data_500, data_600, data_700, data_800, data_1000])
    max_s_vals = np.zeros(len(values))
    p_max = np.zeros(len(values))
    for i in range(0, len(values)):
        max_s_vals[i] = np.max(values[i])
        p_max[i] = np.argmax(values[i])

    return p_max/(len(data_600)), max_s_vals

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



def plot_all_curves():
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(data_100)
    ax.plot(data_200)
    #ax.plot(data_300)
    ax.plot(data_400)
    ax.plot(data_500)
    ax.plot(data_600)
    ax.plot(data_700)
    ax.plot(data_800)
    #ax.plot(data_900)

    plt.figure(2)
    plt.plot(n_sites, np.exp(intercept)*n_sites**slope)
    plt.plot(n_sites, max_s_vals, "x")
    plt.xscale("log")
    plt.yscale("log")
    plt.show()
plot_all_curves()
plot_ny()
