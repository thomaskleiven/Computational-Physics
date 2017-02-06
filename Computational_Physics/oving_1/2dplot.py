import numpy as np
from matplotlib import pyplot as plt


u0 = 1
D = 1
L = 1
start = L/2


def f1(t, x, x0):
    N = 6000
    result = np.zeros(len(x))
    for n in range(N, 0, -1):
        result += u0*np.exp(-(np.pi*n/L)**(2)*D*t)*np.sqrt(2/L)*np.sin(n*np.pi*x/L)*np.sqrt(2/L)*np.sin(n*x0*np.pi/L)
    return result/N




t = np.linspace(0, 0.01, 1000)
x = np.linspace(0, L, 1000)

y1 = f1(1/16000, x, start)
y2 = f1(0.001, x, start)
y3 = f1(0.01, x, start)

print(y1)

plt.plot(x, y1, "bo")
plt.plot(x, y2)
plt.plot(x, y3, "ro")

#plt.xlabel(r'$\tau$')
#plt.ylabel(r'P($\tau$)')
#plt.axis([0.0, 1, 0.0, 2])
plt.grid(True)
plt.show()


#def func(x,t):
#    #return (((u0/np.sqrt(4*np.pi*D*t))*np.exp(-(x-x0)*(x-x0)/(4*D*t))))
#    return u0*(sum( (np.exp(-(n*np.pi/3)**(2)*D*t*np.sqrt(2/3)*np.sin(n*np.pi*x/3)*np.sqrt(2/3)*np.sin(n*np.pi*x0/3))) for n in range(0,1000)))
