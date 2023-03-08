# This code solves the orbit of the Earth around the Sun using the Verlet
# method. The first plot shows the orbit. The second plot shows the potential,
# kinetic, and total energy of the Earth during its orbit. The third plot shows
# just the total energy to illustrate the small oscillation in total energy
# during the Earth's orbit.
#
# d^2r/dt^2 = -GMr/r^3

import numpy as np
import matplotlib.pyplot as plt

G = 6.6738e-11*(60*60)**2 # m3 / kg hr2
M = 1.9891e30 # kg
m = 5.9722e24 # kg

def f(rv,t):
    x,y,vx,vy = rv
    fx = vx
    fy = vy
    R = np.sqrt(x**2+y**2)
    fvx = -G*M*x/(R**3)
    fvy = -G*M*y/(R**3)
    return np.array([fx,fy,fvx,fvy],float)

t_min = 0
t_max = 10*365*24
h = 1
N = int((t_max-t_min)/h)
t_val = np.arange(t_min,t_max,h)

init_x = 1.4710e11 # m
init_y = 0 # m
init_vx = 0 # m / hr
init_vy = 3.0287e4*60*60 # m / hr
rv = np.array([init_x,init_y,init_vx,init_vy],float)
x_val = np.zeros(N)
y_val = np.zeros(N)
vx_val = np.zeros(N)
vy_val = np.zeros(N)
vx_half = init_vx+0.5*h*f(rv,t_val[0])[2]
vy_half = init_vy+0.5*h*f(rv,t_val[0])[3]
v_half = np.array([vx_half,vy_half],float)
pot_NRG = np.zeros(N)
kin_NRG = np.zeros(N)
tot_NRG = np.zeros(N)
for i in range(N):
    x_val[i],y_val[i],vx_val[i],vy_val[i] = rv
    rv[0:1+1] = rv[0:1+1]+h*v_half
    k = h*f(rv,t_val[i])[2:]
    rv[2:] = v_half+0.5*k
    v_half = v_half+k
    R = np.sqrt(rv[0]**2+rv[1]**2)
    pot_NRG[i] = -G*M*m/R
    V = np.sqrt(rv[2]**2+rv[3]**2)
    kin_NRG[i] = 0.5*m*V**2
    tot_NRG[i] = pot_NRG[i]+kin_NRG[i]

print('\na)')
plt.plot(x_val,y_val,label='Earth')
plt.scatter(0,0,color='y',label='Sun')
plt.legend()
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()

print('\nb)')
plt.plot(t_val/24/365,pot_NRG,label='Potential energy')
plt.plot(t_val/24/365,kin_NRG,label='Kinetic energy')
plt.plot(t_val/24/365,tot_NRG,label='Total energy')
plt.legend()
plt.xlabel('Time (years)')
plt.ylabel('Energy')
plt.show()

print('\nc)')
plt.plot(t_val/24/365,tot_NRG,label='Total energy')
plt.legend()
plt.xlabel('Time (years)')
plt.ylabel('Energy')
plt.show()