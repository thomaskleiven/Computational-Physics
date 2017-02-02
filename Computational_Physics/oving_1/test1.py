from math import factorial, cos, e, sqrt
from scipy import *
import numpy as np
from matplotlib import pyplot as plt


def f1(t):
    return 0.5*(1 + sum( (a**(2*n)*cos(2*sqrt(1 + n)*t))/(e**a**2*factorial(n)) for n in range(0,100)))

a=4
t = linspace(0, 35, 1000)
y1 = f1(t)

print(y1)

plt.plot(t, y1)

plt.show()
