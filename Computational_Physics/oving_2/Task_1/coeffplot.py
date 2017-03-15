import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

data = np.loadtxt("coeff.csv", delimiter="\n")
data.astype(int)
#data = np.reshape(data, (3,3))

x = np.zeros(19999)
x[0] = 1/20000
n = int(len(x))
for i in range(1, n):
    x[i] = x[i-1]+1/40000

#print(len(data))
print(x[19998])

plt.plot(x, data)
plt.show()

#plt.imshow(data)
#plt.show()
#plt.plot(data)
#plt.show()
