from matplotlib import pyplot as plt
import numpy as np

solution = np.loadtxt("psi.csv", delimiter=",")
plt.plot(solution)
plt.show()
