from matplotlib import pyplot as plt
import numpy as np
import sys
from os import listdir
from os.path import isfile, join

param = "with"
path = "eigenvectors/"
files = [f for f in listdir(path) if isfile(join(path, f)) and param in f]

try:
    n = int(files[0][27:31])
    print "loaded file with %d"%n + " eigenmodes"
except:
    n = int(files[0][27:30])
    print "loaded file with %d"%n + " eigenmodes"

mode = np.loadtxt("eigenvectors/eigenvector_with_potential_%d.csv"%n, delimiter=",")
mode = np.transpose(mode)

plt.plot(mode[0])
plt.show()
