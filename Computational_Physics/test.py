import numpy as np

a = np.matrix('3 1 1; -1 3 1')

s,u,v = np.linalg.svd(a, full_matrices = True)


print(s)
print(u)
print(v)
