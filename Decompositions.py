# This code contains functions to perform a variety of matrix decompositions:
# - LU decomposition (LUdec)
# - LU decomposition with partial pivoting (LUdecPP)
# - QR decomposition (QRdec)

import numpy as np

def LUdec(A):
    N = len(A[:,0])
    L = np.zeros((N,N))
    U = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            if i <= j:
                U[i,j] = A[i,j]
                for k in range(i): U[i,j] -= U[k,j]*L[i,k]
            if i >= j:
                L[i,j] = A[i,j]
                for k in range(j): L[i,j] -= U[k,j]*L[i,k]
                L[i,j] /= U[j,j]
    return L,U

def LUdecPP(A0):
    A = np.copy(A0)
    N = len(A[:,0])
    swaps = np.zeros(N)
    L = np.zeros((N,N))
    U = np.zeros((N,N))
    for i in range(N):
        row_max = np.where(np.abs(A[i:,i]) == np.max(np.abs(A[i:,i])))[0][0]+i
        A[[i,row_max]] = A[[row_max,i]]
        swaps[i] = row_max
        for j in range(N):
            if i <= j:
                U[i,j] = A[i,j]
                for k in range(i): U[i,j] -= U[k,j]*L[i,k]
            if i >= j:
                L[i,j] = A[i,j]
                for k in range(j): L[i,j] -= U[k,j]*L[i,k]
                L[i,j] /= U[j,j]
    return L,U,swaps.astype(int)

def QRdec(A0):
    N = len(A0[0,:])
    A = np.copy(A0)
    Q = np.zeros((N,N))
    R = np.zeros((N,N))
    for i in range(N):
        u = A[:,i]
        for j in range(i):
            R[j,i] = np.dot(Q[:,j],A[:,i])
            u -= np.dot(R[j,i],Q[:,j])
        R[i,i] = np.linalg.norm(u)
        Q[:,i] = u/R[i,i]
    return Q,R