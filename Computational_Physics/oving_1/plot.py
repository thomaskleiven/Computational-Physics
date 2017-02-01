import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt("implicitEuler.csv", delimiter=",")

x = np.linspace(0.0, 1.0, len(data[:,0]))
plt.plot(x,data)
plt.xlim([0.4990,0.5010])
#plt.ylim([0,100])
plt.show()
