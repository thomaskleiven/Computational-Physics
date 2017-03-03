from matplotlib import pyplot as plt
import numpy as np

data = np.loadtxt("ConservativeLaxWendroff.csv", delimiter=",")
plt.plot(data)
plt.show()
