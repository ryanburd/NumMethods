# This code solves two variable ordinary differential equations using the
# fourth-order Runge-Kutta method. Specifically, the Lotka-Volterra equations
# for the predator-prey interactions between biological species are solved,
# plotting the population of the predator (fox) and prey (rabbit) vs time.
#
# dx/dt = alpha*x - beta*xy, dy/dt = gamma*xy - delta*y

import numpy as np
import matplotlib.pyplot as plt

alpha = 1
beta = 0.5
gamma = 0.5
delta = 2

def f(r,t):
    x,y = r
    fx = alpha*x-beta*x*y
    fy = gamma*x*y-delta*y
    return np.array([fx,fy],float)

t_min = 0
t_max = 30
N = 1000
h = (t_max-t_min)/N
t_val = np.arange(t_min,t_max,h)

init_x = 2
init_y = 2
r = np.array([init_x,init_y],float)
x_val = np.zeros(N)
y_val = np.zeros(N)
for i in range(N):
    x_val[i],y_val[i] = r
    k1 = h*f(r,t_val[i])
    k2 = h*f(r+0.5*k1,t_val[i]+0.5*h)
    k3 = h*f(r+0.5*k2,t_val[i]+0.5*h)
    k4 = h*f(r+k3,t_val[i]+h)
    r += (k1+2*k2+2*k3+k4)/6
plt.plot(t_val,x_val,label='rabbits')
plt.plot(t_val,y_val,label='foxes')
plt.xlabel('Time')
plt.ylabel('Population, in thousands')
plt.legend(bbox_to_anchor=[1.0,1.0])
plt.show()