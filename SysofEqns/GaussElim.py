# This code performs Gaussian elimination with partial pivoting on a linear
# system of equations. The matrix A and vector v must be provided as numpy
# arrays of type float. To obtain the solution from the output, use BackSub in
# Substitutions.py

import numpy as np
from Substitutions import BackSub

def GaussElimPP(A0,v):
    A = np.copy(A0)
    N = len(v)
    for j in range(N):
        row_max = np.where(np.abs(A[j:,j]) == np.max(np.abs(A[j:,j])))[0][0]+j
        A[[j,row_max]] = A[[row_max,j]]
        temp = v[j]
        v[j] = v[row_max]
        v[row_max] = temp
        div = A[j,j]
        A[j,:] /= div
        v[j] /= div
        for i in range(j+1,N):
            mult = A[i,j]
            A[i,:] -= mult*A[j,:]
            v[i] -= mult*v[j]
    return A,v

if __name__ == "__main__":
    A = np.array([[2,  1,  4,  1],
                [3,  4, -1, -1],
                [1, -4,  1,  5],
                [2, -2,  1,  3]], float)
    v = np.array([-4, 3, 9, 7], float)
    A, v = GaussElimPP(A,v)
    b = BackSub(A,v)
    print(b)