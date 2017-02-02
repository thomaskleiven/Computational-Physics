import numpy as np
from matplotlib import pyplot as plt
import pylab as py

data = np.loadtxt("explicitEuler.csv", delimiter=",")

x = np.linspace(0, 1.0, len(data[:,0]))
plt.plot(x,data)
plt.xlim([0.499,0.501])
#plt.ylim([0,100])
plt.show()
