# This code calculates the volume of a 10D unit hypersphere using Monte Carlo
# integration. The volume is printed.

import numpy as np

N = int(1e7)
inside = 0

for i in range(N):
    x = np.random.random(10)*2-1
    sum_x2 = 0
    for j in range(len(x)): sum_x2 += x[j]**2
    if sum_x2 <= 1: inside += 1
integral = 2**len(x)*inside/N

print('Volume = %.3f'%integral)