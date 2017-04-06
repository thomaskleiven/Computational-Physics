from matplotlib import pyplot as plt
import numpy as np

mode = np.loadtxt("euler_scheme.csv", delimiter=",")
plt.plot(mode)
plt.show()
