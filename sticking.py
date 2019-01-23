import sys, math, argparse
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
parser = argparse.ArgumentParser()
def InpParams(parser):
    parser.add_argument("t", help='temperature', default='300', nargs='?')
    args = parser.parse_args()
    if args.t:
        temperature = args.t
    else:
        temperature = '300'
    return int(temperature)

Temp = InpParams(parser)

if Temp == 80:
    exponent = 1.5
elif Temp == 190:
    exponent = 1.0
elif Temp == 300:
    exponent = 0.5
else:
    exponent = 2.0

fname = '/home/becker/lammps/DelayStickingProb.dat'
f = open(fname, 'r')

Arr = np.zeros([3,100])

ctr = 0
for i,line in enumerate(f):
    ang, temp, nrg, stick, sigma = line.split()
    if int(temp) == Temp:
        Arr[:,ctr] = float(nrg)*math.cos(float(ang)*math.pi/180.)**exponent, float(stick), float(sigma)
        ctr += 1

xlbl = r"E $\cos^n(\vartheta$) / meV"
# plt.axis([0,100,0,1])
plt.xlabel(xlbl)
plt.ylabel(r"Sticking Probability $R_{st}$")
# plt.title("T = " + str(Temp) + " K, n = " + str(exponent))
plt.errorbar(Arr[0,:ctr], Arr[1,:ctr], Arr[2,:ctr], ls='None', marker='x', color='black')
plt.tight_layout()
plt.savefig(str(Temp)+ "_" + str(int(exponent*10)) + "InitSticking.pdf")
plt.show()
