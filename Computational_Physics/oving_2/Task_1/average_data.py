import numpy as np
from matplotlib import pyplot as plt
import sys
from os import listdir
from os.path import isfile, join

params = ["avg100", "avg200", "avg300", "avg400", "avg500", "avg600", "avg700", "avg800", "avg900", "avg1000", "avg1500", "avg2000", "p100", "p200", "p300", "p400", "p500", "p600", "p700", "p800", "p900", "p1000", "p1500", "p2000", "chi100", "chi200", "chi300", "chi400", "chi500", "chi600", "chi700", "chi800", "chi900", "chi1000", "chi1500", "chi2000"]

path = sys.argv[1]
#param = sys.argv[2]

i=0
for param in params:
    files = [f for f in listdir(path) if isfile(join(path, f)) and param in f]
    if(len(files) == 0):
        continue

    data = np.zeros(len(np.loadtxt((path+"/"+files[0]), delimiter="\n")[40:-40]))

    for filename in files:
        data += np.loadtxt((path+"/"+filename), delimiter="\n")[40:-40]

    data /= len(files)

    np.savetxt((param+"_averaged.csv"), data, delimiter="\n")
    print "File saved as %s.csv"%param
