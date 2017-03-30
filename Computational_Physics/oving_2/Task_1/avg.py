import numpy as np
from matplotlib import pyplot as plt
import sys

filename = sys.argv[1]

data = np.loadtxt('results/'+filename+'.csv', delimiter='\t')

x = np.linspace(0.0,1.0,len(data))

plt.plot(x, data)
plt.grid()
plt.show()
