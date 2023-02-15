# This code solves Laplace's equation for a 2D square box with voltage V = 1
# volt applied to the top and 0 elsewhere. The Jacobi method is used with speed
# up from the combined overrelaxation/Gauss-Seidel method to solve the partial
# differential equation. Several values of the overrelaxation paramter omega
# are tested to determine the value that requires the least iterations. The
# voltage over of the square box is plotted.
#
# d^2phi/dx^2 + d^2phi/dy^2 + d^2phi/dz^2 = 0

import numpy as np
import matplotlib.pyplot as plt

M = 100
a = 1 #cm
V = 1.0 #V

omega = np.arange(0.8,1,0.025)
counts = np.zeros(len(omega))
for w in range(len(omega)):
    phi = np.zeros((M+1,M+1),float)
    phi[0,:] = V
    eps = 1e-6
    largest_delta = 1e6
    counter = 0
    while largest_delta > eps:
        counter += 1
        largest_delta = 0
        for x in range(M+1):
            for y in range(M+1):
                if x==0 or x==M or y==0 or y==M: pass
                else:
                    delta = phi[x,y]
                    phi[x,y] = (1+omega[w])/4*(phi[x+1,y]+phi[x-1,y]+phi[x,y+1]+phi[x,y-1])-omega[w]*phi[x,y]
                    delta -= phi[x,y]
                    if abs(delta) > abs(largest_delta):
                        largest_delta = abs(delta)
    counts[w] = counter
    print('Omega = %.3f, # of iterations = %d'%(omega[w],counter))
plt.imshow(phi)
plt.xlabel('cm')
plt.ylabel('cm')
plt.colorbar()
plt.show()