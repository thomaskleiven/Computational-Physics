import numpy as np
from matplotlib import pyplot as plt


crank = np.loadtxt("crankNicolson_dx.csv", delimiter=",")
delX = 1.0/(len(crank)+1)

L=0.5
dPlus = 1.0
dMinus = 0.1

t = 0.07
aPlus = 2/(1+ np.sqrt(dMinus/dPlus))
aMinus = aPlus*np.sqrt(dMinus/dPlus)

#uPlus   = aPlus/(np.sqrt(np.pi*4*dPlus*t))*np.exp(-(x*x)/(4*dPlus*t))
#uMinus = aMinus/(np.sqrt(np.pi*4*dMinus*t))*np.exp(-(x*x)/(4*dPMinus*t))

def fPlus(t, x):
    return aPlus/(np.sqrt(np.pi*4*dPlus*t))*np.exp(-(x*x)/(4*dPlus*t))

def fMinus(t, x):
    return aMinus/(np.sqrt(np.pi*4*dMinus*t))*np.exp(-(x*x)/(4*dMinus*t))

xPlus = np.linspace(0, L, 100)
xMinus = np.linspace(-L,0, 100)
x = np.linspace(-L+delX,L-delX, 401)

y1 = fPlus(0.01, xPlus)
y2 = fMinus(0.01, xMinus)

#Plot numerical
plt.figure(2)
plt.plot(x, crank)
plt.plot(x, crank,'o', mfc='none', markevery = 10)
plt.grid(True)


#Plot analytical
plt.plot(xPlus, y1)
plt.plot(xMinus, y2)
plt.grid(True)
plt.show()
