import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

Lx = 1
Ly = 1

#Plot evolution matrix
for i in range(1,100,3):
    matrix = np.loadtxt('data/matrix_%d.csv'%i, delimiter=",")
    plt.imshow(matrix, cmap = "coolwarm")
    plt.show()


#Plot eigenmode
mode = 1
data = np.loadtxt('data/eigenvector_%d.csv'%mode, delimiter=",")

plt.figure(1)
plt.imshow(data, cmap = "coolwarm")
plt.savefig('figures/mode%d.svg'%mode)


#Plot eigenvalues
numEigenvalues = np.loadtxt('data/eigenvalues.csv', delimiter=",")
plt.figure(2)
plt.yscale('log')
plt.plot(numEigenvalues, "+")

#normalize = numEigenvalues.shape[0]

def eigen_val(kx, ky):
    return ((kx)**(2)/Lx**(2) + (ky)**(2)/(Ly)**(2))*np.pi**(2)

eigenvalues = np.zeros(len(numEigenvalues))

eigenvalues[0] = eigen_val(1,1)
eigenvalues[1] = eigen_val(2,1)
eigenvalues[2] = eigen_val(1,2)
eigenvalues[3] = eigen_val(2,2)
eigenvalues[4] = eigen_val(3,1)
eigenvalues[5] = eigen_val(1,3)
eigenvalues[6] = eigen_val(2,3)
eigenvalues[7] = eigen_val(3,2)
eigenvalues[8] = eigen_val(4,1)
eigenvalues[9] = eigen_val(1,4)
eigenvalues[10] = eigen_val(3,3)
eigenvalues[11] = eigen_val(2,4)
eigenvalues[12] = eigen_val(4,2)
eigenvalues[13] = eigen_val(3,4)
eigenvalues[14] = eigen_val(4,3)
eigenvalues[15] = eigen_val(5,1)
eigenvalues[16] = eigen_val(1,5)
eigenvalues[17] = eigen_val(5,2)
eigenvalues[18] = eigen_val(2,5)
eigenvalues[19] = eigen_val(4,4)
eigenvalues[20] = eigen_val(5,3)
eigenvalues[21] = eigen_val(3,5)





print (eigenvalues)
plt.plot(eigenvalues, "*", color = 'red')



plt.show()
