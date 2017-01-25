import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import csv
import numpy as np
import math as mp
from scipy import stats

x = []
y = []

with open('data.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        y.append(np.float(row[0]))


x = np.arange(len(y))

#plt.hist(y, bins = 10000)
#plt.xlim([0,2000])
#plt.ylim([0,30000])
#plt.yscale('log')
#plt.figure(2)
#plt.xscale('log')
#plt.show()


#Plot cumultative distributin
values, base = np.histogram(y, bins = 10000, normed=1)
plt.yscale('log')
plt.xscale('log')
plt.plot(base[:-1], values, c='green', ls="dashdot")
plt.xlim([0,100000])

print (len(base), len(values))

y_min = 500
y_max = 10000
base_min = np.min(base)
base_max = np.max(base)


values = values[base < y_max]
base = base[base < y_max]
values = values[base > y_min]
base = base[base > y_min]

slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(base), np.log(values))

x = np.linspace(base_min, base_max, 100)
fitted_line = np.exp(intercept)*x**slope
print(slope)
plt.plot(x, fitted_line)





plt.show()




#Calculate alpha
#list_of_logs = []
#for i in range(0,270):
#    list_of_logs.append(np.float(-np.log(cumulative[i]) - np.float(np.log(base[i]))))
#
#alpha = sum(list_of_logs)/len(list_of_logs)
#print (alpha)
