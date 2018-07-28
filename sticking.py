import sys, math
import matplotlib.pyplot as plt
import numpy as np
Temp = 190

if Temp == 80:
    exponent = 1.5
elif Temp == 190:
    exponent = 1.0
elif Temp == 300:
    exponent = 0.5
else:
    exponent = 2.0

fname = 'testdata.txt'
f = open(fname, 'r')

Arr = np.zeros([3,100])

ctr = 0
for i,line in enumerate(f):
    ang, temp, nrg, stick, sigma = line.split()
    if int(temp) == Temp:
        print('Hallo')
        Arr[:,ctr] = float(nrg)*math.cos(float(ang)*math.pi/180.)**exponent, float(stick), float(sigma)
        ctr += 1


plt.errorbar(Arr[0,:ctr], Arr[1,:ctr], Arr[2,:ctr])
plt.show()
