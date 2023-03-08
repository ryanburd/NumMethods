# This code solves the wave equation for a piano string being struck at a
# position d from the end using the forward-time centered-space (FTCS) method,
# plotting an animation of the wave in time until the solution of the partial
# differential equation becomes numerically unstable.
#
# d^2phi/dx^2 - 1/v^2*d^2phi/dt^2 = 0
#
# As FTCS equations:
# phi(x,t+h) = phi(x,t) + h*psi(x,t)
# psi(x,t+h) = psi(x,t) + h*v^2/a^2*[phi(x+a,t) + phi(x-a,t) - 2phi(x,t)]

import numpy as np
from vpython import *

C = 1 # m / s
L = 1 # m
d = 0.1 # m
sigma = 0.3 # m

def psi_init(x):
    return C*x*(L-x)/L**2*np.exp(-(x-d)**2/(2*sigma**2))

v = 100 # m / s
h = 1e-5 # s
a = 0.05 # m
total_time = 0.7 # s
time_points = int(total_time/h)

x_val = np.arange(0,L+a/2,a)
phi = np.zeros(len(x_val))
psi = psi_init(x_val)
phi_prev = phi
vbstr = curve(origin=vector(-0.5,0,0))
for i in range(len(x_val)):
    vbstr.append(vector(x_val[i],phi[i],0))
for t in range(time_points):
    for i in range(1,len(x_val)-1):
        phi[i] += h*psi[i]
        psi[i] += h*(v**2)/(a**2)*(phi_prev[i+1]+phi_prev[i-1]-2*phi_prev[i])
        vbstr.modify(i,y=1e3*phi[i])
    print(h*t)