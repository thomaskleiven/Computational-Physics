import numpy as np
from matplotlib import pyplot as plt
import sys
from os import listdir
from os.path import isfile, join

#params = ["avg100_", "avg200_", "avg300_", "avg400_", "avg500_", "avg600_", "avg700_", "avg800_", "avg900_", "avg1000_", "avg1500_", "avg2000_"]
#params = ["p100_", "p200_", "p300_", "p400_", "p500_", "p600_", "p700_", "p800_", "p900_", "p1000_", "p1500_", "p2000_"]
params = ["chi100_", "chi200_", "chi300_", "chi400_", "chi500_", "chi600_", "chi700_", "chi800_", "chi900_", "chi1000_", "chi1500_", "chi2000_"]

path = sys.argv[1]
#param = sys.argv[2]

i=0
for param in params:
    files = 0
    files = [f for f in listdir(path) if isfile(join(path, f)) and param in f]
    if(len(files) == 0):
        continue

    data = np.zeros(len(np.loadtxt((path+"/"+files[0]), delimiter="\n")[40:-40]))

    for filename in files:
        data += np.loadtxt((path+"/"+filename), delimiter="\n")[40:-40]

    data /= len(files)
    print "File averaged: ", len(files)

    np.savetxt((param+"_averaged.csv"), data, delimiter="\n")
    print "File saved as %s.csv"%param
