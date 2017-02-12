import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats

mode = 4
data = np.loadtxt("data/eigenvector_%d_50.csv"%mode, delimiter = ",")

ny = data.shape[0]
nx = data.shape[1]

data *= np.sqrt(nx*ny)
data = -data

n = 2 #Mode 2,2

dx = 1.0/(nx+1)     #Excluded boundry
dy = 1.0/(ny+1)

kx = n*np.pi
ky = n*np.pi

x = np.linspace(dx,1.0-dx, nx)
y = np.linspace(dy, 1.0-dy, ny)

Lx = 1
Ly = 1

nDisc = [6, 12, 25, 50, 100]
errorNorm = 2



def errorScaling(mode):
    errors = np.zeros(5)
    for i in range(0,5):
        fname = "data/eigenvector_%d_%d.csv"%(mode, nDisc[i])
        num = np.loadtxt(fname, delimiter = ",")

        ny = num.shape[0]
        nx = num.shape[1]

        num *= np.sqrt(nx*ny)

        dx = 1.0/(nx+1)     #Excluded boundry
        dy = 1.0/(ny+1)

        kx = 2*np.pi
        ky = 2*np.pi

        x = np.linspace(dx,1.0-dx, nx)
        y = np.linspace(dy, 1.0-dy, ny)
        X,Y = np.meshgrid(x,y)
        amplitude = eigen_func(X,Y, kx, ky)

        """plt.figure(1)
        plt.imshow(num, cmap = "coolwarm")
        plt.colorbar()

        plt.figure(2)
        X,Y = np.meshgrid(x,y)
        plt.imshow(amplitude, cmap = "coolwarm")
        plt.colorbar()
        plt.show()"""

        errors[i] = error(num, amplitude)

    return errors



def eigen_func(x, y, kx, ky):
    return 2/(np.sqrt(Lx/Ly))*np.sin(kx*x/Lx)*np.sin(ky*y/Ly)

def eigen_val(kx, ky):
    return ((kx)**(2)/Lx**(2) + (ky)**(2)/(Ly)**(2))*np.pi**(2)


def error(numerical, analytical):
    normalization = numerical.shape[0] * numerical.shape[1]
    return (np.sum(np.abs(np.abs(numerical) - np.abs(analytical))**errorNorm))**(1.0/errorNorm)/normalization

slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(nDisc), np.log(errorScaling(4)))

print(slope)

plt.plot(nDisc, errorScaling(4), color = "black")
plt.plot(nDisc, errorScaling(4), 'o', mfc='none')
plt.xlabel("Number of discretication points")
plt.ylabel("Error")
plt.xscale('log', base = 10)
plt.yscale('log', base = 10)
plt.show()


X,Y = np.meshgrid(x,y)
amplitude = eigen_func(X,Y, kx ,ky)


errMat = np.abs(data - amplitude)

plt.figure(1)
plt.imshow(errMat, cmap = "nipy_spectral", norm = mpl.colors.LogNorm())
plt.colorbar()

plt.figure(2)
plt.imshow(data, cmap = "coolwarm")
plt.colorbar()

plt.figure(3)
plt.imshow(amplitude, cmap = "coolwarm")
plt.colorbar()

plt.show()
