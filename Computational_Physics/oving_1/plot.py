import matplotlib.pyplot as plt
import csv
import numpy as np

x = []
y = []

with open('data.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        y.append(np.float(row[0]))

x = range(0,len(y))

'''
values, base = np.histogram(y, bins = 40)

cumulative = np.cumsum(values)

# plot the cumulative function
plt.plot(base[:-1], cumulative, c='blue')
#plot the survival function
plt.plot(base[:-1], len(y)-cumulative, c='green')

plt.show()'''


plt.axhline(linewidth=1, color= 'r')

plt.plot(x,y, label='x(t) = ∑s_n')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Interesting Graph\nCheck it out')
plt.legend()
plt.show()
