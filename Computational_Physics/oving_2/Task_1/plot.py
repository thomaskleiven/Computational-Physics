import numpy as np
from matplotlib import pyplot as plt
import sys
from os import listdir
from os.path import isfile, join

path = sys.argv[1]
param = sys.argv[2]

files = [f for f in listdir(path) if isfile(join(path, f)) and param in f]

data = np.zeros(len(np.loadtxt((path+"/"+files[0]), delimiter="\n")[40:-40]))

print len(files)

for filename in files:
    data += np.loadtxt((path+"/"+filename), delimiter="\n")[40:-40]

data /= len(files)
np.savetxt(("averaged"+files[0]+".csv"), data, delimiter="\n")
print "Averaged file saved"

x = np.linspace(0.0,1.0,len(data))

one = np.loadtxt((path+"/"+files[1]), delimiter="\n")[40:-40]

plt.figure(2)
plt.plot(x, one)

plt.figure(1)
plt.plot(x, data)
plt.show()
