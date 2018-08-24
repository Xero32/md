
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
kB = 1.3806503e-23
e0 = 1.60217662e-19


def expo(x,t,u):
    return t * np.exp((u * 1e-3 * e0) / (x * kB))

def expcorr(x,t,u,a):
    return t * np.exp((u * 1e-3 * e0) / (x * kB) - a)

tmax = 1.1
tmid = 1.2
tmin = 1.3
umax = 86.0
umid = 80.0
umin = 74.0

def test(x,t,u):
    return t * np.exp((u * 1e-3 * e0) / (x * kB))

fname = "ResidenceTimeNorm.dat"
f = open(fname, 'r')
data = np.full((3,30), 999999.9)

for i, line in enumerate(f):
    if len(line.split()) == 0:
        break
    ang, temp, effTemp, energy, tR, sigma = line.split()

    data[0,i] = float(effTemp)
    data[1,i] = float(tR)
    data[2,i] = float(sigma) * 1.2

xErr = [5.0 for i in range(len(data[0]))]
TArr = np.arange(50,1000,1)
p0 = (1.5,80)
p1 = (1.5,80,0)
plt.axis([100,300,0,1500])
eopt, ecov = curve_fit(expo, data[0], data[1], p0=p0, bounds=([0,70],[1.5,100]), sigma=data[2])
copt, ccov = curve_fit(expcorr, data[0], data[1], p0=p1, bounds=([0,70,0],[5,120,10],), sigma=data[2])
plt.errorbar(data[0], data[1], data[2], xerr=xErr, fmt='x', color='black')
#plt.plot(TArr, expo(TArr, *eopt), 'b-', label=r'fit: $%5.3f \cdot exp\left({\frac{%5.3f meV}{k_BT}}\right)$' % tuple(eopt))
#plt.plot(TArr, expcorr(TArr, *copt), 'r-', label=r'fit: $%5.3f \cdot exp\left({\frac{%5.3f meV}{k_BT}}-%5.3f\right)$' % tuple(copt))
plt.plot(TArr, test(TArr, tmax, umax), 'r--', label=r'Max: $t^{R}(T) = %5.1f \cdot exp\left({\frac{%5.1f meV}{k_BT}}\right)$' %(tmax, umax))
plt.plot(TArr, test(TArr, tmid, umid), 'g-', label=r'Mid: $t^{R}(T) = %5.1f \cdot exp\left({\frac{%5.1f meV}{k_BT}}\right)$' %(tmid, umid))
plt.plot(TArr, test(TArr, tmin, umin), 'b-.', label=r'Min: $t^{R}(T) = %5.1f \cdot exp\left({\frac{%5.1f meV}{k_BT}}\right)$' %(tmin, umin))
plt.xlabel('T / K')
plt.ylabel(r'$t^{R}$ / ps')
plt.legend(prop={'size': 12})
plt.tight_layout()
plt.savefig("ResidenceTimeFit_xerr.pdf")
plt.show()
