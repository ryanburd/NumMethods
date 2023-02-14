# This code performs the QR algorithm to obtain the eigenvalues and
# eigenvectors as a matrix A. A must be provided as a numpy array. The function
# outputs the vector d, which contain the eigenvalues, and the matrix V, which
# contains the eigenvectors as columns. The jth eigenvalue corresponds to the
# jth eigenvector.

import numpy as np
from Decompositions import QRdec

def eigQR(A):
    D = np.copy(A)
    N = len(A[0,:])
    V = np.identity(N)
    eps = 1e-6
    test = True
    while test:
        Q,R = QRdec(D)
        D = np.dot(R,Q)
        V = np.dot(V,Q)
        test = False
        for i in range(N):
            for j in range(N):
                if i!=j and abs(D[i,j]) > eps:
                    test = True
                    break
    d = np.zeros(N)
    for i in range(N):
        d[i] = D[i,i]
    return d, V

A = np.array([[1, 2],
              [2, 1]], float)
d, V = eigQR(A)
print(d)
print(V)