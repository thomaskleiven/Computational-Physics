from matplotlib import pyplot as plt
import numpy as np
from scipy import optimize, fmin

v0 = 1000

p = np.linspace(0.0, 1.0, 10000)

def eigvecs(x):
    k = np.sqrt(x)
    kappa = np.sqrt(v0-x)
    first = np.exp(kappa/3)
    second = (kappa*np.sin(k/3) + kappa*np.cos(k/3))**2
    third = np.exp(-kappa/3)
    fourth = (kappa*np.sin(k/3) - kappa*np.cos(k/3))**2

    return second - (1.0/first)*third*fourth

#Find minimum
xopt = optimize.fmin(eigvecs, 271)
print "Optimal x-value: ", xopt[0]

#Find root
res = optimize.brentq(eigvecs, 272.030991745, 273)
print "Root: ", res
