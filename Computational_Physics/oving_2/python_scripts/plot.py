import numpy as np
from matplotlib import pyplot as plt
import sys
from os import listdir
from os.path import isfile, join

path = sys.argv[1]


data = np.loadtxt((path), delimiter="\n")[40:-40]
x = np.linspace(0.0,1.0,len(data))

plt.figure(1)
plt.plot(x, data)
plt.show()
