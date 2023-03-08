# This code solves a second-order ODE as two coupled first-order ODEs,
# specifically that of an (an)harmonic oscillator. The first plot shows a
# harmonic oscillator with 2 different amplitudes. The second plot shows an
# anharmonic oscillator with 3 different amplitudes. The third (fourth) plot
# shows the phase space of the harmonic (anharmonic) oscillator with 2
# different amplitudes. The fourth plot shows the phase space of the van der
# Pol oscillator with 3 different values of mu in:
#
# Harmonic osc: d^2x/dt^2 = -w^2x
# Anharmonic osc: d^2x/dt^2 = -w^2x^3
# van der Pol osc: d^2x/dt^2 - mu(1-x^2)dx/dt + w^2x = 0

import numpy as np
import matplotlib.pyplot as plt

omega = 1

def har_osc(r,t):
    x,y = r
    fx = y
    fy = -omega**2*x
    return np.array([fx,fy],float)

def anhar_osc(r,t):
    x,y = r
    fx = y
    fy = -omega**2*x**3
    return np.array([fx,fy],float)

def vanderPol(r,t):
    x,y = r
    fx = y
    fy = mu*(1-x**2)*y-omega**2*x
    return np.array([fx,fy],float)

t_min = 0
t_max = 50
N = 1000
h = (t_max-t_min)/N
t_val = np.arange(t_min,t_max,h)

def pendulum(f,init_x):
    init_y = 0
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
    return x_val,y_val

plt.plot(t_val,pendulum(har_osc,1)[0],label='a) Amplitude = 1')
plt.plot(t_val,pendulum(har_osc,2)[0],label='b) Amplitude = 2')
plt.title('Harmonic oscillator pendulum')
plt.xlabel('t')
plt.ylabel('x')
plt.legend()
plt.show()

plt.plot(t_val,pendulum(anhar_osc,1)[0],label='c) Amplitude = 1')
plt.plot(t_val,pendulum(anhar_osc,2)[0],label='c) Amplitude = 2')
plt.plot(t_val,pendulum(anhar_osc,0.5)[0],label='c) Amplitude = 0.5')
plt.title('Anharmonic oscillator pendulum')
plt.xlabel('t')
plt.ylabel('x')
plt.legend()
plt.show()

plt.plot(pendulum(har_osc,1)[0],pendulum(har_osc,1)[1],label='d) Amplitude = 1')
plt.plot(pendulum(har_osc,2)[0],pendulum(har_osc,2)[1],label='d) Amplitude = 2')
plt.title('Harmonic pendulum phase space')
plt.xlabel('x')
plt.ylabel('dx/dt')
plt.legend()
plt.show()

plt.plot(pendulum(anhar_osc,1)[0],pendulum(anhar_osc,1)[1],label='d) Amplitude = 1')
plt.plot(pendulum(anhar_osc,2)[0],pendulum(anhar_osc,2)[1],label='d) Amplitude = 2')
plt.title('Anharmonic pendulum phase space')
plt.xlabel('x')
plt.ylabel('dx/dt')
plt.legend()
plt.show()

t_min = 0
t_max = 20
N = 10000
h = (t_max-t_min)/N
t_val = np.arange(t_min,t_max,h)

mu = 1
plt.plot(pendulum(vanderPol,1)[0],pendulum(vanderPol,1)[1],label='e) $\mu$ = 1')
mu = 2
plt.plot(pendulum(vanderPol,1)[0],pendulum(vanderPol,1)[1],label='e) $\mu$ = 2')
mu = 4
plt.plot(pendulum(vanderPol,1)[0],pendulum(vanderPol,1)[1],label='e) $\mu$ = 4')
plt.title('van der Pol pendulum phase space')
plt.xlabel('x')
plt.ylabel('dx/dt')
plt.legend()
plt.show()