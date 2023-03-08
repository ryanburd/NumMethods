# This code solves the orbit of the Earth (first plot) and Pluto (second plot)
# using the Bulirsch-Stoer method.
#
# d^2x/dt^2 = -GMx/r^3, d^2y/dt^2 = -GMy/r^3

import numpy as np
import matplotlib.pyplot as plt

G = 6.6738e-11*(60*60*24)**2 # m3 / kg days2
M = 1.9891e30 # kg

def f(r,t,x0,y0):
    x,v = r
    fx = v
    R = np.sqrt(x0[0]**2+y0[0]**2)
    fv = -G*M*x/(R**3)
    return np.array([fx,fv],float)

def orbit(init_x,init_y,init_vx,init_vy):
    x0 = np.array([init_x,init_vx],float)
    y0 = np.array([init_y,init_vy],float)
    x_val = np.zeros(N)
    y_val = np.zeros(N)
    
    delta = 1000/7 # m / day
    eps = H*delta
    for i in range(N):
        x_val[i],y_val[i] = x0[0],y0[0]
        n = 1
        Rx_prev = np.zeros((n,2))
        Ry_prev = np.zeros((n,2))
        Rx_error = np.array([H*delta*10,0])
        Ry_error = np.array([H*delta*10,0])
        while abs(Rx_error[0]) > eps or abs(Ry_error[0]) > eps:
            h = H/n
            Rx_cur = np.zeros((n,2))
            Ry_cur = np.zeros((n,2))
            x1,y1 = x0+0.5*h*f(x0,t_val[i],x0,y0),y0+0.5*h*f(y0,t_val[i],x0,y0)
            x2,y2 = x0+h*f(x1,t_val[i],x0,y0),y0+h*f(y1,t_val[i],x0,y0)
            for j in range(1,n):
                x1 += h*f(x2,t_val[i],x0,y0)
                y1 += h*f(y2,t_val[i],x0,y0)
                x2 += h*f(x1,t_val[i],x0,y0)
                y2 += h*f(y1,t_val[i],x0,y0)
            Rx_cur[0],Ry_cur[0] = 0.5*(x2+x1+0.5*h*f(x2,t_val[i],x0,y0)),0.5*(y2+y1+0.5*h*f(y2,t_val[i],x0,y0))
            for j in range(1,n):
                Rx_error = (Rx_cur[j-1]-Rx_prev[j-1])/((n/(n-1))**(2*j)-1)
                Rx_cur[j] = Rx_cur[j-1]+Rx_error
                Ry_error = (Ry_cur[j-1]-Ry_prev[j-1])/((n/(n-1))**(2*j)-1)
                Ry_cur[j] = Ry_cur[j-1]+Ry_error
            Rx_error,Ry_error = Rx_error,Ry_error
            Rx_prev,Ry_prev = Rx_cur,Ry_cur
            n += 1
        x0,y0 = Rx_cur[-1],Ry_cur[-1]
    return x_val,y_val

print('\na)')
t_min = 0 # days
t_max = 365*100 # days
H = 7 # days
N = int((t_max-t_min)/H)
t_val = np.arange(t_min,t_max,H)
init_x = 1.4710e11 # m
init_y = 0 # m
init_vx = 0 # m / day
init_vy = 3.0287e4*60*60*24 # m / day
x_val,y_val = orbit(init_x,init_y,init_vx,init_vy)
plt.plot(x_val,y_val,label='Earth')
plt.scatter(0,0,label='Sun',color = 'y')
plt.legend()
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()

print('\nb)')
plt.plot(x_val,y_val,label='Earth')
plt.scatter(0,0,label='Sun',color = 'y')
t_min = 0 # days
t_max = 365*1000 # days
H = 365 # days
N = int((t_max-t_min)/H)
t_val = np.arange(t_min,t_max,H)
init_x = 4.4368e12 # m
init_y = 0 # m
init_vx = 0 # m / day
init_vy = 6.1218e3*60*60*24 # m / day
x_val,y_val = orbit(init_x,init_y,init_vx,init_vy)
plt.plot(x_val,y_val,label='Pluto')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.legend()
plt.show()