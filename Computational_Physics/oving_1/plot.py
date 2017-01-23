import matplotlib.pyplot as plt
import csv
import numpy as np

x = []
y = []

with open('data.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        y.append((row))

x = range(0,len(y))

plt.axhline(linewidth=1, color= 'r')

plt.plot(x,y, label='Loaded from file!')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Interesting Graph\nCheck it out')
plt.legend()
plt.show()
