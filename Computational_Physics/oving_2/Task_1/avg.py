import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt('avg.csv', delimiter='\t')

x = np.linspace(0.0,1.0,len(data))

plt.plot(x, data)
plt.show()
