# This code uses Monte Carlo integration with importance sampling to calculate
# an integral:
#
# I = 1/N*sum[f(x_i)/w(x_i)*int0_1(w(x)dx)]
# f(x) = x^(-1/2)/(e^x-1)
# w(x) = x^(-1/2)
#
# Histograms of the random numbers and the sample points x are plotted, showing
# that the sample points match the expected probability distribution also
# plotted. The integral solution is then printed, expected to be ~0.84 for this
# example.

import numpy as np
import matplotlib.pyplot as plt

print('\nProblem 3: Exercise 10.8')

def x(z):
    return z**2

def f_over_w(xi):
    return (np.exp(xi)+1)**(-1)

N = int(1e6)
z_val = np.zeros(N)
x_val = np.zeros(N)
sumN = 0
for i in range(N):
    z_val[i] = np.random.random()
    x_val[i] = x(z_val[i])
    sumN += f_over_w(x_val[i])

fig,ax = plt.subplots(1,2)
ax[0].hist(z_val,bins=150)
ax[0].set_xlabel('z')
ax[1].hist(x_val,bins=150,label='x(z)=z$^2$')
ax[1].set_xlabel('x')
ax[1].legend()

def p(t):
    return 1/(2*t**0.5)
t_val = np.arange(0,1,1e-3)
sy = ax[1].twinx()
sy.plot(t_val,p(t_val),color='orange',label=r'p(x)=$\frac{1}{2\sqrt{x}}$')
sy.legend(bbox_to_anchor=(1,0.9))

print('\na)')
plt.show()

integral = sumN*2/N
print('\nb) Integral = %.3f'%integral)