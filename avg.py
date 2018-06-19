import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x,a,b,c):
    return a / (x**b) + c

def inv(x,a,b,c):
    return (a / (x-c)) ** (1./b)

def expo(x,a,b,c):
    return a*np.exp(-x*b)+c

def maximum(a,b,c=1e8):
    if a > b:
        if a > c:
            return c
        else:
            return a
    else:
        if b > c:
            return c
        else:
            return b



def minimum(a,b):
    if a < b:
        return a
    else:
        return b


def sort(Arr):
    NewArr = np.zeros([len(Arr)])
    upper = max(Arr)
    lower = 0
    val = 0
    for i in range(len(NewArr)):

        for j in range(0, len(Arr)):
            if lower < Arr[j] and Arr[j] < upper:
                mini = minimum(v)

        NewArr[i] = mini



lines = 100
data = np.zeros([8, lines])
fname = '/home/becker/lammps/ResTimeStd.dat'

f = open(fname, 'r')

Temps = set()

for i, line in enumerate(f):
    try:
        data[:,i] = line.split()
        Temps.add(data[1,i])
    except:
        continue

f.close()
nd = np.zeros([3, len(Temps)-1])
temps = list()
for i in range(len(Temps)):
    temps.append(Temps.pop())

temps.remove(80)

for t in range(len(temps)):
    ctr = 0
    avg = 0
    std = 0
    for i in range(0, lines):
        if data[1,i] == temps[t]:
            avg += data[3,i]
            std = maximum(std, data[4,i], 10000)
            ctr += 1
    avg /= ctr
    nd[0,t] = temps[t]
    nd[1,t] = avg
    nd[2,t] = std
'''
fnamenew = '/home/becker/lammps/ResTimeAvg.dat'
fn = open(fnamenew, 'w')
for i in range(0, len(temps)):
    fn.write("%d %f %f\n" %(nd[0,i], nd[1,i], nd[2,i]))

fn.close()
'''
#sys.exit()
p0 = (5000, 0.01, 80)
X = np.arange(0,2100,1)
poptb,pcovb= curve_fit(func, nd[1], nd[0], bounds=([0,0.4999999,0],[2500,0.5,100]))
popt, pcov = curve_fit(func, nd[1], nd[0])
eopt, ecov = curve_fit(expo, nd[1], nd[0], p0=p0, sigma=nd[2], bounds=([500, 0.001, 0],[1000000000, 1, 100]))
plt.title('Residence Time')
plt.errorbar(nd[1], nd[0], xerr=nd[2], fmt='x')
plt.axis([1,2000,100,500])
plt.plot(X, func(X, *popt), 'r-', label='fit: %5.3f / t^(%5.3f) + %5.3f' % tuple(popt))
plt.plot(X, func(X, *poptb), 'k-', label='bounded fit: %5.3f / t^(%5.3f) + %5.3f' % tuple(poptb))
plt.plot(X, expo(X, *eopt), 'g-', label='fit: %5.3f * exp(-x * %5.3f) + %5.3f' % tuple(eopt))
plt.xlabel('tR / ps')
plt.ylabel('T / K')
plt.legend()
#plt.savefig('/home/becker/lammps/tmpplot/ResTimeInv.pdf', format='pdf')
plt.show()

######  same procedure, but inverted ######

TArr = np.arange(0,500,1)
poptb,pcovb= curve_fit(inv, nd[0], nd[1], bounds=([0,0.49999999,0],[2500,0.5,100]), sigma=nd[2])
popt, pcov = curve_fit(inv, nd[0], nd[1])
eopt, ecov = curve_fit(expo, nd[0], nd[1], p0=p0, bounds=([0, 0.0001, 0],[1e10, 1, 80]), sigma=nd[2])
plt.title('Residence Time')
plt.errorbar(nd[0], nd[1], yerr=nd[2], fmt='x')
plt.axis([0,350,0,6000])
plt.plot(TArr, inv(TArr, *poptb), 'k-', label='bounded fit: (%5.3f / [T - %5.3f])^%1.0f' %(poptb[0], poptb[2], 1/poptb[1]))
plt.plot(TArr, expo(TArr, *eopt), 'g-', label='fit: %5.3f * exp(-x * %5.3f) + %5.3f' % tuple(eopt))
plt.xlabel('T / K')
plt.ylabel('tR / ps')
plt.legend()
plt.tight_layout()
#plt.savefig('/home/becker/lammps/tmpplot/ResTimeExpo.pdf', format='pdf')
plt.show()
