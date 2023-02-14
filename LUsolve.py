# This code solves linear systems of equations using LU decomposition with
# (LUPPsolve) or without (LUsolve) partial pivoting. The decompositions,
# backward, and forward subsitions are performed in the imported functions. The
# matrix A and vector v must be provided as numpy arrays.

import numpy as np
from Decompositions import LUdec,LUdecPP
from Substitutions import BackSub,ForSub

def LUsolve(A,v):
    L,U = LUdec(A)
    y = ForSub(L,v)
    x = BackSub(U,y)
    return x

def LUPPsolve(A,v):
    L,U,swaps = LUdecPP(A)
    for i in range(len(v)):
        temp = v[i]
        v[i] = v[swaps[i]]
        v[swaps[i]] = temp
    y = ForSub(L,v)
    x = BackSub(U,y)
    return x

A = np.array([[2,  1,  4,  1],
              [3,  4, -1, -1],
              [1, -4,  1,  5],
              [2, -2,  1,  3]], float)
v = np.array([-4, 3, 9, 7], float)
b = LUsolve(A,v)
bPP = LUPPsolve(A,v)
print(b)
print(bPP)