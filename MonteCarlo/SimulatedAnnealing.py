# This code finds the global minimum of a function using simulated annealing.
# Two examples functions are used:
#
# f(x) = x^2 - cos(4pi*x)
# f(x) = cos(x) + cos(sqrt(2)x) + cos(sqrt(3)x) for 0<x<50
#
# For each example, the function and the values of x at each iteration vs time
# are plotted. The convergence of x at the global minimum can be easily seen.

import numpy as np
import matplotlib.pyplot as plt

# Example 1

def fa(x):
    return x**2-np.cos(4*np.pi*x)

Tmax = 10
Tmin = 1e-3
tau = 1e4

x = 2
x_val = [2]
F = fa(x)

t = 0
t_val = [0]
T = Tmax
while T>Tmin:
    t += 1
    t_val.append(t)
    T = Tmax*np.exp(-t/tau)
    z = np.random.random()
    r = np.sqrt(-2*np.log(1-z))
    theta = np.random.random()*2*np.pi
    xnew = x+r*np.cos(theta)
    Fnew = fa(xnew)
    dF = Fnew-F
    if np.random.random()<np.exp(-dF/T):
        x = xnew
        F = Fnew
    x_val.append(x)

fig,ax = plt.subplots(1,2)
x = np.arange(-2,2,0.01)
ax[0].plot(x,fa(x))
ax[0].set_xlabel('x')
ax[0].set_ylabel('f(x)')
ax[1].plot(t_val,x_val,".",markersize=1)
ax[1].set_xlabel('Time')
ax[1].set_ylabel('x')
plt.show()

# Example 2

def fb(x):
    return np.cos(x)+np.cos(np.sqrt(2)*x)+np.cos(np.sqrt(3)*x)

Tmax = 10
Tmin = 1e-3
tau = 1e6

x = 2
x_val = [2]
F = fb(x)

t = 0
t_val = [0]
T = Tmax
while T>Tmin:
    t += 1
    t_val.append(t)
    T = Tmax*np.exp(-t/tau)
    z = np.random.random()
    r = np.sqrt(-2*np.log(1-z))
    theta = np.random.random()*2*np.pi
    xnew = x+r*np.cos(theta)
    while xnew<=0 or xnew>=50:
        z = np.random.random()
        r = np.sqrt(-2*np.log(1-z))
        theta = np.random.random()*2*np.pi
        xnew = x+r*np.cos(theta)
    Fnew = fb(xnew)
    dF = Fnew-F
    if np.random.random()<np.exp(-dF/T):
        x = xnew
        F = Fnew
    x_val.append(x)

fig,ax = plt.subplots(1,2)
x = np.arange(0,50,0.01)
ax[0].plot(x,fb(x))
ax[0].set_xlabel('x')
ax[0].set_ylabel('f(x)')
ax[1].plot(t_val,x_val,".",markersize=1)
ax[1].set_xlabel('Time')
ax[1].set_ylabel('x')
plt.show()