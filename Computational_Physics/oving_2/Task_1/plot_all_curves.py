import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

data_70 = np.loadtxt('res_hon/p70.csv', delimiter='\n')[10:-10]
data_300 = np.loadtxt('res_hon/p300.csv', delimiter='\n')[10:-10]
data_400 = np.loadtxt('res_hon/p400.csv', delimiter='\n')[10:-10]
data_500 = np.loadtxt('res_hon/p500.csv', delimiter='\n')[10:-10]
data_600 = np.loadtxt('res_hon/p600.csv', delimiter='\n')[10:-10]
data_700 = np.loadtxt('res_hon/p700.csv', delimiter='\n')[10:-10]
data_800 = np.loadtxt('res_hon/p800.csv', delimiter='\n')[10:-10]
data_900 = np.loadtxt('res_hon/p900.csv', delimiter='\n')[10:-10]
data_1000= np.loadtxt('res_hon/p1000.csv', delimiter='\n')[10:-10]
data_5000= np.loadtxt('res_hon/p1500.csv', delimiter='\n')[10:-10]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
x = np.linspace(0.0, 1.0, len(data_900))

ax.plot(x, data_300)
ax.plot(x, data_400)
ax.plot(x, data_500)
ax.plot(x, data_600)
ax.plot(x, data_700)
ax.plot(x, data_800)
ax.plot(x, data_900)
ax.plot(x, data_1000)
ax.plot(x, data_5000)

plt.show()
