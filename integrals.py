# This code contains functions to solve integrals numerically with a variety of
# methods:
# - Trapezoidal rule (adap_trap)
# - Simpson's rule (adap_Simpson)
# - Romberg
# - Gaussian quadrature (GaussQuad)
#
# Use of GaussQuad requires gaussxw.

import numpy as np
from gaussxw import gaussxw

def adap_trap(f,a,h,N,old_I):
    new_sum = 0
    for k in range(1,N,2): new_sum += f(a+k*h)
    new_I = 0.5*old_I+h*new_sum
    error = 1/3*(new_I-old_I)
    return new_I,error

def adap_Simpson(f,a,h,N,old_I,old_S,old_T):
    new_S = old_S+old_T
    sum_odd = 0
    for k in range(1,N,2): sum_odd += f(a+k*h)
    new_T = 2/3*sum_odd
    new_I = h*(new_S+2*new_T)
    error = 1/15*(new_I-old_I)
    return new_I,new_S,new_T,error

def Romberg(f,a,h,N,old_I,old_R):
    new_sum = 0
    for k in range(1,N,2): new_sum += f(a+k*h)
    new_I = 0.5*old_I+h*new_sum
    new_R = np.zeros(len(old_R)+1)
    new_R[0] = new_I
    for m in range(1,len(old_R)+1):
        error = 1/(4**m-1)*(new_R[m-1]-old_R[m-1])
        new_R[m] = new_R[m-1]+error
    return new_I,new_R,error

def GaussQuad(f,a,b,N):
    x,w = gaussxw(N)
    xr = 0.5*(b-a)*x+0.5*(b+a)
    wr = 0.5*(b-a)*w
    integral = 0
    for k in range(N): integral += wr[k]*f(xr[k])
    return integral