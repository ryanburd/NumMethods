# This code provides an example of using the Gaussian quadrature method to
# calculate the period of an anharmonic oscillator as a function of the
# amplitude, a. The potential is V(x) = x^4.

import numpy as np
from gaussxw import gaussxw
import matplotlib.pyplot as plt

print('Problem 3: Exercise 5.10')

def f(x,a):
    return 1/np.sqrt(a**4-x**4)

def period(a):
    m = 1
    N = 20
    x,w = gaussxw(N)
    xr = 0.5*a*x+0.5*a
    wr = 0.5*a*w
    integral = 0
    for k in range(N): integral += wr[k]*f(xr[k],a)
    return np.sqrt(8*m)*integral

inc = 0.01
a_values = np.arange(0,2+inc/2,inc)
period_values = np.zeros(len(a_values))
for a in range(len(a_values)):
    period_values[a] = period(a_values[a])
plt.plot(a_values,period_values)
plt.xlabel('Amplitude, a')
plt.ylabel('Period, T')
plt.show()