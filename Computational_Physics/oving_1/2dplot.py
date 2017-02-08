import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

#data = np.loadtxt("explicitEuler.csv", delimiter=",")

explicitEuler = False
implicitEuler = False
crankNicolson = False


impEuler5 = np.loadtxt("crankNicolson.csv", delimiter=",")

#data = np.loadtxt("explicitEuler.csv", delimiter=",")
#data = np.loadtxt("implicitEuler.csv", delimiter=",")


u0 = 1
D = 1
L = 1
delX = 1.0/(len(impEuler5)+1)
start = (int(len(impEuler5)/2)+1)*delX


if(explicitEuler):
    x = np.linspace(0, L, len(impEuler5))
else:
    x = np.linspace(delX, 1-delX, len(impEuler5))


def f1(t, x, x0, reflective):
    N = 90000
    result = np.zeros(len(x))
    for n in range(N, 0, -1):
        if(reflective == False):
            result += u0*np.exp(-(np.pi*n/L)**(2)*D*t)*np.sqrt(2/L)*np.sin(n*np.pi*x/L)*np.sqrt(2/L)*np.sin(n*x0*np.pi/L)
        else:
            if(n == 1):
                result += np.sqrt(1/L)
                continue
            result += u0*np.exp(-(np.pi*n/L)**(2)*D*t)*np.sqrt(2/L)*np.cos(n*np.pi*x/L)*np.sqrt(2/L)*np.cos(n*x0*np.pi/L)

    return result

def error(exact, numerical):
    return np.sum(np.abs(exact - numerical))/len(numerical)



y1 = f1(0.002, x, start, True)


if(explicitEuler):
    euler23 = np.loadtxt("explicitEuler_33.csv", delimiter=",")
    euler46 = np.loadtxt("explicitEuler_3.csv", delimiter=",")
    euler92 = np.loadtxt("explicitEuler_0.csv", delimiter=",")
    #euler184 = np.loadtxt("explicitEuler_3.csv", delimiter=",")
    deltaT = 1/np.array([3E4, 3E5, 3E6])
    errorsExplicit = [error(y1, euler23), error(y1, euler46), error(y1, euler92)]
    #linregress
    slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(deltaT), np.log(errorsExplicit))

    #Plot Errors
    plt.figure(1)
    #plt.plot(deltaT, errorsExplicit)
    plt.plot(deltaT, np.exp(intercept)*deltaT**slope)
    plt.plot(deltaT, errorsExplicit,'o', mfc='none')
    plt.xscale('log', basex = 2)
    plt.yscale('log', basey = 10)
    plt.grid(True)
elif(implicitEuler):
    impEuler10 = np.loadtxt("implicitEuler_100.csv", delimiter=",")
    impEuler21 = np.loadtxt("implicitEuler_50.csv", delimiter=",")
    impEuler43 = np.loadtxt("implicitEuler_25.csv", delimiter=",")
    impEuler5 = np.loadtxt("implicitEuler_500.csv", delimiter=",")
    #Timesteps
    deltaT = 1/np.array([2E3, 1E4, 4E4])
    errorsImplicit = [error(y1, impEuler5), error(y1, impEuler10),error(y1, impEuler43)]
    #linregress
    slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(deltaT), np.log(errorsImplicit))
    #Plot Errors
    plt.figure(1)
    #plt.plot(deltaT, errorsImplicit)
    plt.plot(deltaT, np.exp(intercept)*deltaT**slope)
    plt.plot(deltaT, errorsImplicit,'o', mfc='none')
    plt.xscale('log', basex = 2)
    plt.yscale('log', basey = 10)
    plt.grid(True)
elif(crankNicolson):
    #crank500 = np.loadtxt("crankNicolson_2.csv", delimiter=",")
    crank100 = np.loadtxt("crankNicolson_20.csv", delimiter=",")
    crank50 = np.loadtxt("crankNicolson_200.csv", delimiter=",")
    crank10 = np.loadtxt("crankNicolson_2000.csv", delimiter=",")
    #Timesteps
    deltaT = 1/np.array([5E2, 5E3, 5E4])
    errorsCrank = [error(y1, crank10), error(y1, crank50), error(y1, crank100)]
    #linregress
    slope, intercept, rvalue, pvalue, stdr = stats.linregress(np.log(deltaT), np.log(errorsCrank))
    #Plot Errors
    plt.figure(1)
    #plt.plot(deltaT, errorsImplicit)
    plt.plot(deltaT, np.exp(intercept)*deltaT**slope)
    plt.plot(deltaT, errorsCrank,'o', mfc='none')
    plt.xscale('log', basex = 2)
    plt.yscale('log', basey = 10)
    plt.grid(True)



#Plot the functions
plt.figure(2)
plt.plot(x, impEuler5)
plt.plot(x, impEuler5,'o', mfc='none', markevery = 15)
plt.plot(x, y1)
#plt.plot(x, y2)
#plt.plot(x, y3)
plt.grid(True)
plt.show()
