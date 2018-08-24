######################################################################
# Program to analyze equilibration of Surface for MD simulations
# usage: python3 sth.py [temperature] []

import matplotlib
import sys
import math
import matplotlib.pyplot as plt
from array import array
import numpy as np
import smooth
import histogram as h
import plot
import matplotlib
matplotlib.rcParams.update({'font.size': 16})

checkname = sys.argv[2]
temp = sys.argv[1]
T = int(temp)

j=0
kB = 1.3806503e-23
e0 = 1.60217662e-19 #(in Coulomb)
amu = 1.66054e-27 #(amu to kg)
#m = 197.0*amu
m = 195.0*amu #Pt
e0inv = 1 / e0
Q = np.zeros([6,10000])
V = np.zeros([3,10000])
VV = np.zeros([3,10000])
E = np.zeros([5,10000])

######################
# distributions
######################

def boltz1(m,t,E):
    return np.exp(-E/(kB*t))/( 2. * np.sqrt(np.pi * kB * t * E))

def boltz2(m,t,E):
    return (1./(kB*t))*np.exp(-E/(kB*t))

def boltz3(m,t,E):
    return 2. * np.sqrt(E/np.pi) * (1/(kB*t))**(3./2.) * np.exp(-E/(kB*t))

def boltz3v(m,t,v):
    return (m/(2.*np.pi*kB*t))**(3./2.)*4.*np.pi*v**2 * np.exp(-m*v**2/(2.*kB*t))
    #return np.sqrt(E/np.pi) * np.exp(-E/(kB*t))

def pKE(vx,vy): #[J]
    return 0.5 * m * (vx**2 + vy**2)

def nKE(vz): #[J]
    return 0.5 * m * vz**2

def Gauss(t, mu, h):
    g = np.exp( (-1.0) * (t-mu)*(t-mu) / (2*h*h) )
    return g

def GaussSmoothing(N, t, f, xarr):
    h = 5
    fG = 0
    GaussSum = 0
    for j in range (0,N):
        GaussSum += Gauss(t, xarr[j], h)
    for j in range (0,N):
        fG += Gauss(t, xarr[j], h) / GaussSum * f[j]
    return fG # double

########################

fl_q = open(checkname,'r')

for i, line in enumerate(fl_q):
    if i < 1:
        continue
    else:
        if math.fabs(float(line.split()[6])) <= 10e-7 and math.fabs(float(line.split()[7])) <= 10e-7 and math.fabs(float(line.split()[8])) <= 10e-7:
            continue
        else:
            iden, temp, pe, Q[0,j], Q[1,j], Q[2,j], Q[3,j], Q[4,j], Q[5,j] = line.split()
            j += 1

VV[0,:j] = Q[5,:j] * 100
VV[1,:j] = np.sqrt((100*Q[3,:j])**2 + (100*Q[4,:j])**2)
VV[2,:j] = np.sqrt((100*Q[3,:j])**2 + (100*Q[4,:j])**2 + (100*Q[5,:j])**2)
c = min(VV[2,:j])
d = max(VV[2,:j])

#T = 300.0

#T = float(temp)
E[0,:j] = nKE( Q[3,:j] * 100)
E[1,:j] = nKE( Q[4,:j] * 100)
E[2,:j] = nKE( Q[5,:j] * 100)
E[3,:j] = pKE( Q[3,:j] * 100, Q[4,:j] * 100)
E[4,:j] = E[2,:j] + E[3,:j]
a=min(E[3,:j])
b=max(E[4,:j])


e = np.arange(a, b, (b-a)/j)
v = np.arange(c,d,(d-c)/j)
out = np.zeros([5,100])
sze = np.zeros([5,100])
for i in range(0,5):
    out[i,:], sze[i,:]  = h.hist(E[i,:j],100, sze[i,:])


#fl_out = open('/home/becker/lammps/eq/eqxl_ns.txt','w')
#fl_velo= open('/home/becker/lammps/eq/eqvelos0.txt','w')
#for i in range(0,j):
    ##fl_out.write('%g %g %g %g %g\n' %(E[0,i], E[1,i], E[2,i], E[3,i], E[4,i]))
    #fl_velo.write('%g %g %g\n' %(VV[0,i], VV[1,i], VV[2,i]))

# plot v_xyz
plt.plot(v, boltz3v(m, T, v), color='black')
quan = r"$v$"
unit = r"$\frac{m}{s}$"
unitinv = unit + r'$^{-1}$'
XLABEL = quan + " / " + unit
YLABEL = r"P("+quan+") / " + unitinv
plt.xlabel(quan + " / " + unit)
plt.ylabel(r"P("+quan+") / " + unitinv)
plt.tight_layout()
if T < 190:
    upper = 250
    plot.Histogram(0,upper,100, VV[2,:j], [],[], lbl=r'$P\left(v=\sqrt{v_x^2+v_y^2+v_z^2}\right)$', Title=str(T)+' K', edge='pos', acol='red', bcol='red', ccol='red', totcol='red',
    saveflag=1, savename='velocityEq.pdf', xlbl=XLABEL, ylbl=YLABEL)
else:
    plot.Histogram(0,400,100, VV[2,:j], [],[], lbl=r'$P\left(v=\sqrt{v_x^2+v_y^2+v_z^2}\right)$', Title=str(T)+' K', edge='pos', acol='red', bcol='red', ccol='red', totcol='red',
    saveflag=1, savename='velocityEq.pdf', xlbl=XLABEL, ylbl=YLABEL)


# plot z Energy
plt.plot(e, boltz1(m, T, e), color='black')
quan = r"$E_z$"
unit = r"$J$"
unitinv = unit + r'$^{-1}$'
plt.xlabel(quan + " / " + unit)
plt.ylabel(r"P("+quan+") / " + unitinv)
if T < 190:
    upper = 8e-21
    plot.Histogram(0, upper, 100, E[2,:j], [], [], lbl='z Energy / J', Title=str(T) + ' K', edge='pos')
else:
    plot.Histogram(0, 1.5e-20, 100, E[2,:j], [], [], lbl='z Energy / J', Title=str(T) + ' K', edge='pos')

#plot xy Energy
plt.plot(e, boltz2(m, T, e), color='black')
quan = r"$E_{xy}$"
unit = r"$J$"
unitinv = unit + r'$^{-1}$'
plt.xlabel(quan + " / " + unit)
plt.ylabel(r"P("+quan+") / " + unitinv)
if T < 190:
    upper = 1e-20
    plot.Histogram(0, upper, 100, E[3,:j], [], [], lbl='xy Energy / J', Title=str(T) + ' K', edge='pos')
else:
    plot.Histogram(0, 2e-20, 100, E[3,:j], [], [], lbl='xy Energy / J', Title=str(T) + ' K', edge='pos')

#plot xyz Energy
plt.plot(e, boltz3(m, T, e), color='black')
quan = r"$E_{xyz}$"
unit = r"$J$"
unitinv = unit + r'$^{-1}$'
plt.xlabel(quan + " / " + unit)
plt.ylabel(r"P("+quan+") / " + unitinv)
if T < 190:
    upper = 1e-20
    plot.Histogram(0, upper, 100, E[4,:j], [], [], lbl='xyz Energy / J', Title=str(T) + ' K', edge='pos')
else:
    plot.Histogram(0, 3e-20, 100, E[4,:j], [], [], lbl='xyz Energy / J', Title=str(T) + ' K', edge='pos')

#fl_out.close()
fl_q.close()
