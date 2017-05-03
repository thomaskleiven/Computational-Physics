from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import sys
from os import listdir
from os.path import isfile, join

import matplotlib as mpl
mpl.rcParams['svg.fonttype'] = "none"
mpl.rcParams['font.size'] = 18

path = "/home/tkleiven/Documents/Skole/6_semester/Computational_Physics/QM/eigenvectors/"
files = [f for f in listdir(path) if isfile(join(path, f))]
eigen_error = np.zeros(len(files))
steps = np.zeros(len(files))

#Analytical eigenfunction
def eigenfunction(n,x):
    return np.sqrt(2)*np.sin((n+1)*np.pi*x)

#Squared error, numerical vs. analytical
def error(numerical, analytical):
    return np.sum((np.abs(numerical) - np.abs(analytical))**2)

#Errorscaling as a function of discretisation for a given eigenvector
def errorScaling():
    for i in range(0,len(files)):
        eigenvector = np.loadtxt(("/home/tkleiven/Documents/Skole/6_semester/Computational_Physics/QM/eigenvectors/"+files[i]), delimiter=",")
        eigenvector = eigenvector[:,0]
        N=len(eigenvector)
        dx = 1.0/N
        x = np.linspace(0.0+dx,1.0-dx,N)
        steps[i] = N
        analytical = eigenfunction(0,x)
        eigenvector *= np.sum(analytical)/np.sum(eigenvector)
        plt.figure(i)
        plt.plot(x, np.abs(eigenvector))
        plt.plot(x, np.abs(analytical))
        eigen_error[i] = np.sqrt(error(eigenvector, analytical))/N
    return steps, eigen_error


def plotOneEigenvec(eigvec_num):
    filename = "/home/tkleiven/Documents/Skole/6_semester/Computational_Physics/QM/eigenvectors/eigenvector_1000.csv"
    eigenvector = np.loadtxt(filename, delimiter=",")
    eigenvector = np.transpose(eigenvector)
    N = len(eigenvector[eigvec_num])
    dx = 1.0/N
    x = np.linspace(0.0+dx, 1.0-dx, N)
    fig = plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    ax.plot(x, -np.sqrt(N)*eigenvector[eigvec_num],"o", label="\$Numerical$", markevery=30, mfc="none")
    ax.plot(x, eigenfunction(eigvec_num, x), color="black", label="\$Analytical")

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.set_xlabel("\$ x (a.u) $")
    ax.set_ylabel("\$ \psi (a.u)$")
    ax.legend(loc="upper right", frameon=False, labelspacing=0.05)

#Errorscaling for a given eigenvector
steps, eigen_error = errorScaling()

#Linear regression
slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(steps), np.log(eigen_error))
print (slope)

a, V = np.polyfit(np.log(steps), np.log(eigen_error), 1, cov=True)
print a
print V


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.figure(2)
ax.plot(steps, np.exp(intercept)*steps**slope, color="black", label="\$Fit a \cdot x+b$")
ax.plot(steps, eigen_error, "x", label="data")
ax.set_xscale('log')
ax.set_yscale('log')
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")
ax.set_xlabel("\$ Number\ of\ discretisation\ steps $")
ax.set_ylabel("\$ Error $")
ax.legend(loc="upper right", frameon=False, labelspacing=0.05)

plt.show()
