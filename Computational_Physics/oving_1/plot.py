import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import csv
import numpy as np
import math as mp

x = []
y = []

with open('data.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        if(np.float(row[0]) == 100000):
            continue
        y.append(np.float(row[0]))


x = np.arange(len(y))


#Plot cumultative distributin
values, base = np.histogram(y, bins = 10000, normed=1)
plt.yscale('log')
plt.xscale('log')
plt.plot(base[:-1], values, c='green', ls="dashdot")
plt.xlim([0,100000])
plt.show()




#Calculate alpha
#list_of_logs = []
#for i in range(0,270):
#    list_of_logs.append(np.float(-np.log(cumulative[i]) - np.float(np.log(base[i]))))
#
#alpha = sum(list_of_logs)/len(list_of_logs)
#print (alpha)
