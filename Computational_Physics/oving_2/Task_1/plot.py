import numpy as np
from matplotlib import pyplot as plt
import sys

filename = sys.argv[1]

data = np.loadtxt(filename+'.csv', delimiter='\t')
data = data[:-20]

x = np.linspace(0.0,1.0,len(data))

plt.plot(x, data)

plt.show()
