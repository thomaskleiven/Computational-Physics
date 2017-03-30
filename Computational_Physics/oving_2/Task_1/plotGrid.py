import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt("SquareLattice.csv", delimiter=",")

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_axis_bgcolor('black')

for i in range(0,data.shape[0]):
    ax.plot([data[i,0], data[i,2]], [data[i,1], data[i,3]], color='#7FFF00', lw=2)

plt.show()
