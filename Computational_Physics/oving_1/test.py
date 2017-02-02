#from scitools.std import *
from math import factorial, cos, e, sqrt
from scipy import *
import numpy as np
from matplotlib import pyplot as plt

u = 1
L = 35
D = 1

data = np.loadtxt("explicitEuler.csv", delimiter=",")



x = linspace(-1,1,100)
x0 = 0.5

def f1(t):
    return sum(exp(-(n*pi/L)*t)*sqrt(2/L)*sin(n*pi*x/L) for n in range(1,400))

t = linspace(0, 1, 100)
y1 = f1(t)

plt.xlim([-1,1])


plt.plot(t, y1)
plt.show()
