# This code uses Monte Carlo to simulate the decay of (213)Bi to (209)Bi
# through two pathways:
# 98%: (213)Bi --> (209)Pb --> (209)Bi
#  2%: (213)Bi --> (209)Tl --> (209)Pb --> (209)Bi
#
# The number of each atom is plotted vs time

import numpy as np
import matplotlib.pyplot as plt

time = np.arange(0,20000,1) #s

Bi213 = np.zeros(len(time),int)
Bi213[0] = 10000
Bi213_tau = 46*60 #s
pBi213 = 1-2**(-1/Bi213_tau)

Tl209 = np.zeros(len(time),int)
Tl209_tau = 2.2*60 #s
pTl209 = 1-2**(-1/Tl209_tau)

Pb209 = np.zeros(len(time),int)
Pb209_tau = 3.3*60 #s
pPb209 = 1-2**(-1/Pb209_tau)

Bi209 = np.zeros(len(time),int)

for t in range(1,len(time)):
    Pb209_decay = 0
    for p in range(Pb209[t-1]):
        if np.random.random() < pPb209:
            Pb209_decay += 1
    Pb209[t] = Pb209[t-1]-Pb209_decay
    Bi209[t] = Bi209[t-1]+Pb209_decay
    
    Tl209_decay = 0
    for p in range(Tl209[t-1]):
        if np.random.random() < pTl209:
            Tl209_decay += 1
    Tl209[t] = Tl209[t-1]-Tl209_decay
    Pb209[t] += Tl209_decay
    
    Bi213_decay = 0
    Tl209_gain = 0
    Pb209_gain = 0
    for p in range(Bi213[t-1]):
        if np.random.random() < pBi213:
            Bi213_decay += 1
            if np.random.random() < 0.0209: Tl209_gain += 1
            else: Pb209_gain += 1
    Bi213[t] = Bi213[t-1]-Bi213_decay
    Tl209[t] += Tl209_gain
    Pb209[t] += Pb209_gain

plt.plot(time,Bi213,label='$^{213}$Bi')
plt.plot(time,Tl209,label='$^{209}$Tl')
plt.plot(time,Pb209,label='$^{209}$Pb')
plt.plot(time,Bi209,label='$^{209}$Bi')
plt.xlabel('Time (seconds)')
plt.ylabel('Number of atoms')
plt.legend()
plt.show()