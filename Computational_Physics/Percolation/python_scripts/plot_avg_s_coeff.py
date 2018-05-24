import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

import matplotlib as mpl
mpl.rcParams['svg.fonttype'] = "none"
mpl.rcParams['font.size'] = 16
mpl.rcParams['legend.numpoints'] = 1

start = 0
end = 10
p_c = 0.499932975871
#p_c = 0.653619302949
#p_c = 0.346481233244
n_sites = np.array([300,500,700,800,1000,2000])

data_300 = np.loadtxt('square/avg300__averaged.csv', delimiter='\n')[40:-40]
data_500 = np.loadtxt('square/avg500__averaged.csv', delimiter='\n')[40:-40]
data_700 = np.loadtxt('square/avg700__averaged.csv', delimiter='\n')[40:-40]
data_800 = np.loadtxt('square/avg800__averaged.csv', delimiter='\n')[40:-40]
data_1000 = np.loadtxt('square/avg1000__averaged.csv', delimiter='\n')[40:-40]
data_2000 = np.loadtxt('square/avg2000__averaged.csv', delimiter='\n')[40:-40]

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
    print ("Slope of |pmax-pc| as a function of correlation length: ", slope)

    #Fit line
    a, V = np.polyfit(np.log(n_sites[start:end]), np.log(p_max_p_c[start:end]), 1, cov=True)

    #print (a)
    #print (V)

    fig = plt.figure(10)
    ax = fig.add_subplot(1,1,1)
    ax.plot(n_sites, np.exp(intercept)*n_sites**slope, label="\$Fit a \cdot x+b$")
    ax.plot(n_sites, p_max_p_c, 'x', label="\$data$")
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.set_xlabel("\$Number of nodes $")
    ax.set_ylabel("\$Average cluster size$")
    ax.legend(loc="upper right", frameon=False, labelspacing=0.05)

    plt.show()

#Fit line to max(<s>)
slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(n_sites[start:end]), np.log(max_s_vals[start:end]))
print "Slope of max(<s>) as a function of correlation length: ", slope

#Fit line
a, V = np.polyfit(np.log(n_sites[start:end]), np.log(max_s_vals[start:end]), 1, cov=True)
#print a, V


fig = plt.figure(9)
ax = fig.add_subplot(1,1,1)
ax.plot(n_sites, np.exp(intercept)*n_sites**slope, label="\$Fit a \cdot x+b$")
ax.plot(n_sites, max_s_vals, 'x', label="\$data$")
ax.set_xscale("log")
ax.set_yscale("log")

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")
ax.set_xlabel("\$Number of nodes $")
ax.set_ylabel("\$log(max(<s>))$")
ax.legend(loc="lower right", frameon=False, labelspacing=0.05)


#a, V = np.polyfit(np.log(n_sites[start:end]), np.log(max_s_vals[start:end]), 1, cov=True)

plot_ny()
