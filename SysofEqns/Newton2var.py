# This code uses Newton's method to solve simultaneous nonlinear equations with
# 2 variables. One root is printed.
#
# Example used:
# x^2 + y^2 - 3 = 0
# xy - 1 = 0

import numpy as np
import matplotlib.pyplot as plt
from GaussElim import GaussElimPP
from Substitutions import BackSub

def f1(x,y):
    return x**2+y**2-3

def f2(x,y):
    return x*y-1

def df1dx(x,y):
    return 2*x

def df1dy(x,y):
    return 2*y

def df2dx(x,y):
    return y

def df2dy(x,y):
    return x

eps = 1e-10
x = 0.5
y = 1.5
J = np.zeros((2,2))
F = np.zeros((2,1))
DX = np.array([x,y],float)

while abs(DX[0]) > eps or abs(DX[1]) > eps:
    J[0,0] = df1dx(x,y)
    J[0,1] = df1dy(x,y)
    J[1,0] = df2dx(x,y)
    J[1,1] = df2dy(x,y)
    F[0,0] = f1(x,y)
    F[1,0] = f2(x,y)
    JR,FR = GaussElimPP(J,F)
    DX = BackSub(JR,FR)
    x -= DX[0]
    y -= DX[1]

print('One root is (x = %.3f, y = %.3f)'%(x,y))