# This code performs backward (BackSub) and forward substitution (ForSub) on a
# provided matrix A and vector v, returning the solution x.

import numpy as np

def BackSub(A,v):
    N = len(v)
    x = np.zeros(N)
    for i in range(1,N+1):
        x[-i] = v[-i]
        for j in range(i-1,0,-1): x[-i] -= A[-i,-j]*x[-j]
        x[-i] /= A[-i,-i]
    return x

def ForSub(A,v):
    N = len(v)
    x = np.zeros(N)
    for i in range(N):
        x[i] = v[i]
        for j in range(i): x[i] -= A[i,j]*x[j]
        x[i] /= A[i,i]
    return x