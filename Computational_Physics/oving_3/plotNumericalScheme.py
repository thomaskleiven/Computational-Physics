from matplotlib import pyplot as plt
import numpy as np

mode = np.loadtxt("crankNicolsonScheme.csv", delimiter=",")
mode = np.transpose(mode)
plt.plot(mode[1])
plt.show()
