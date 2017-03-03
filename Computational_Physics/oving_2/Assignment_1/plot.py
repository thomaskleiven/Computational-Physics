import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

data = np.loadtxt("Sites.csv", delimiter="\n")
data.astype(int)
#data = np.reshape(data, (3,3))

N = np.sqrt(len(data))
positions = np.array([None]*len(data))

x = np.linspace(0, 1)
y = np.linspace(0, 1)


def makeLattice():
    for i in range(0,len(data)):
        x = int(data[i]/N)
        y = int(data[i]%N)
        positions[i] = (x,y)
makeLattice()

X,Y = np.meshgrid(x,y)

print(X)

plt.imshow(X,Y)
plt.show()

#plt.imshow(data)
#plt.show()
#plt.plot(data)
#plt.show()
