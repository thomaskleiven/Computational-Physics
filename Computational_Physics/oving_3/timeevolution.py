from matplotlib import pyplot as plt
import numpy as np
import matplotlib.animation as animation
from pylab import *
import subprocess

probabilitifunctions = np.loadtxt("crankNicolsonScheme.csv", delimiter=",")
probabilitifunctions = np.transpose(probabilitifunctions)
#make folder
subprocess.Popen(["mkdir", "evolution"], stdout=subprocess.PIPE).communicate()[0]

i = 0
for probability in probabilitifunctions:
    plt.ylim([1E-7,10])
    plt.plot(probability)
    plt.savefig("evolution/img"+"%.d"%i+".png")
    plt.close()
    i+=1


#make movie
subprocess.Popen(["ffmpeg", "-r", "13", "-i", "evolution/img%d.png",  "movie.mp4"], stdout=subprocess.PIPE).communicate()[0]
#delete folder
subprocess.Popen(["rm", "-r", "evolution"], stdout=subprocess.PIPE).communicate()[0]
