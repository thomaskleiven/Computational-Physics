import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt("crankNicolson.csv", delimiter=",")
#data = np.loadtxt("explicitEuler.csv", delimiter=",")
#data = np.loadtxt("implicitEuler.csv", delimiter=",")


u0 = 1
D = 1
L = 1
start = L/2
delX = 1.0/len(data)


def f1(t, x, x0, reflective):
    N = 60000
    result = np.zeros(len(data))
    for n in range(N, 0, -1):
        if(reflective == False):
            result += u0*np.exp(-(np.pi*n/L)**(2)*D*t)*np.sqrt(2/L)*np.sin(n*np.pi*x/L)*np.sqrt(2/L)*np.sin(n*x0*np.pi/L)
        else:
            if(n == 0):
                results += np.sqrt(1/L)
                continue
            result += u0*np.exp(-(np.pi*n/L)**(2)*D*t)*np.sqrt(2/L)*np.cos(n*np.pi*x/L)*np.sqrt(2/L)*np.cos(n*x0*np.pi/L)

    return result

def error(exact, numerical):
    return np.sum(np.abs(exact - numerical))

x = np.linspace(delX, 1-delX, len(data[0:]))
y1 = f1(0.07, x, start, False)
y2 = f1(0.01, x, start, False)
y3 = f1(0.04, x, start, False)

#error = error(y1, data)

print(y1)

#Plot the functions
#plt.plot(x, data)
plt.plot(x,data,'o', mfc='none', markevery = 15)
plt.plot(x, y1)
plt.plot(x, y2)
plt.plot(x, y3)
plt.grid(True)
plt.show()
