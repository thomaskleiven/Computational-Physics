from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import sys
from os import listdir
from os.path import isfile, join

path = "eigenvectors/"
files = [f for f in listdir(path) if isfile(join(path, f))]
eigen_error = np.zeros(len(files))
steps = np.zeros(len(files))


def eigenfunction(n,x):
    return np.sqrt(2)*np.sin((n+1)*np.pi*x)

def error(numerical, analytical):
    return np.sum(np.abs(np.abs(numerical) - np.abs(analytical)))**2

def errorScaling(eigvec_num):
    for i in range(0,len(files)):
        eigenvector = np.loadtxt(("eigenvectors/"+files[i]), delimiter=",")
        N=len(eigenvector[0])
        dx = 1.0/N
        x = np.linspace(0.0+dx,1.0-dx,N)
        steps[i] = N
        eigen_error[i] = error(np.sqrt(N)*eigenvector[eigvec_num], eigenfunction(eigvec_num,x))
    return steps, eigen_error

def plotOneEigenvec(eigvec_num):
    filename = path+files[3]
    eigenvector = np.loadtxt(filename, delimiter=",")
    N = len(eigenvector[eigvec_num])
    dx = 1.0/N
    x = np.linspace(0.0+dx, 1.0-dx, N)
    plt.plot(np.sqrt(N) * eigenvector[eigvec_num])
    plt.plot(eigenfunction(eigvec_num,x))
    plt.figure(1)

plotOneEigenvec(0)

steps, eigen_error = errorScaling(3)
slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(steps), np.log(eigen_error))
print (slope)

analytical = np.loadtxt("analyticalEigenvector.csv", delimiter=",");
plt.figure(4)
plt.plot(analytical)

plt.figure(2)
plt.plot(steps, np.exp(intercept)*steps**slope)
plt.plot(steps, eigen_error, "x")
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Number of disc step")
plt.ylabel("Error")
plt.show()
