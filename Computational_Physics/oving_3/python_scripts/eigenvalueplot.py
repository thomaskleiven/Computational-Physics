from matplotlib import pyplot as plt
import numpy as np

data = np.loadtxt("eigenvalues.csv", delimiter=",")
x = np.zeros(100)

for i in range(0,100):
    x[i] = (i*np.pi)**2

plt.plot(x, data, "x")
plt.show()
