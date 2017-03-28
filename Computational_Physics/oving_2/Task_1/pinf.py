import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt('p.csv', delimiter='\t')

x = np.linspace(0.4,0.6,len(data))

plt.plot(x, data)
plt.show()
