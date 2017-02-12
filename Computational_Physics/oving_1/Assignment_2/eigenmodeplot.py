import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


#Plot eigenmodes
mode = 100
data = np.loadtxt('data/eigenvector_%d.csv'%mode, delimiter=",")

plt.figure(1)
plt.imshow(data, cmap = "coolwarm")
plt.savefig('figures/mode%d.svg'%mode)


#Plot eigenvalues

data = np.loadtxt('data/eigenvalues.csv', delimiter=",")
plt.figure(2)
plt.yscale('log')
plt.plot(data, "+")



plt.show()
