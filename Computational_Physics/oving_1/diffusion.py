import numpy
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D

csvfile = 'crankNicolson.csv'#input("Velg cvsfil:\n")
data = numpy.loadtxt(csvfile, delimiter=",")

figure = pyplot.figure();
ax = figure.add_subplot(111, projection='3d')

#x = numpy.linspace(0, len(data[0,:]), len(data[:,0]))
X, Y = numpy.meshgrid(range( len(data[0,:] )), range( len(data[:,0])  ))

ax.set_zlim3d(0, .2)
ax.plot_surface(X, Y, data)
#pyplot.plot(x,data)
#pyplot.xlim([0.499,0.501])
#plt.ylim([0,100])
pyplot.show()
