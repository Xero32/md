import numpy as np
import plot
import matplotlib.pyplot as plt

pi = 3.14159265359
TempArr = [190, 240, 300]
AngArr = [0,30,45,60]
BncArr = [1,2,4,6,8,10,12,14]
size = 35

def Assertion(ang, temp, bnc):
    assert(bnc==1 or bnc==2 or bnc==4 or bnc==6 or bnc==8 or bnc==10 or bnc==12 or bnc==14)
    assert(temp==80 or temp==120 or temp==160 or temp==190 or temp==240 or temp==300)
    assert(ang==0 or ang==30 or ang==45 or ang==60)

def TempCompare(Angle, Bounce, Theta):
    XArr = []
    fl = open('/home/becker/lammps/111/AngDistr.dat', 'r')
    ctr = 0
    for i, line in enumerate(fl):
        a, t, e, nb, x, theta = line.split()
        if i < 35:
            XArr.append(float(x))
        # look at what effect does the surface temperature have on the ang distr
        if float(a) == Angle and float(nb) == Bounce:
            Theta[:,ctr] = float(a),float(t),float(e),int(nb),float(x),float(theta)
            ctr += 1

    for num, temp in enumerate(TempArr):
        PlotArr = []
        for i in range(len(Theta[0,:])):
            if Theta[-1,i] > 0:
                if Theta[1,i] == TempArr[num]:
                    PlotArr.append(Theta[5,i])

        plt.plot(XArr, PlotArr[:size], label=str(temp)+' K')

    fl.close()
    locs = [0, pi/12, pi/6, pi/4, pi/3, 5*pi/12, pi/2]
    labels = ['0', r'$\frac{\pi}{12}$', r'$\frac{\pi}{6}$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{3}$', r'$\frac{5\pi}{12}$', r'$\frac{\pi}{2}$']
    plt.xticks(locs, labels)
    plot.SetupPlot('Angle: '+str(Theta[0,0])+', Bounce: '+str(Theta[3,0]), r'$\theta$ / rad', r'p($\theta)$ / rad$^{-1}$')



def AngleCompare(Temp, Bounce, Theta):
    XArr = []
    fl = open('/home/becker/lammps/111/AngDistr.dat', 'r')
    ctr = 0
    for i, line in enumerate(fl):
        a, t, e, nb, x, theta = line.split()
        if i < 35:
            XArr.append(float(x))
        # look at what effect does the incident angle have on the ang distr
        if float(t) == Temp and float(nb) == Bounce:
            Theta[:,ctr] = float(a),float(t),float(e),int(nb),float(x),float(theta)
            ctr += 1


    for num, ang in enumerate(AngArr):
        PlotArr = []
        for i in range(len(Theta[0,:])):
            if Theta[-1,i] > 0:
                if Theta[0,i] == AngArr[num]:
                    PlotArr.append(Theta[5,i])

        plt.plot(XArr, PlotArr[:size], label=str(ang)+' deg')

    fl.close()
    locs = [0, pi/12, pi/6, pi/4, pi/3, 5*pi/12, pi/2]
    labels = ['0', r'$\frac{\pi}{12}$', r'$\frac{\pi}{6}$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{3}$', r'$\frac{5\pi}{12}$', r'$\frac{\pi}{2}$']
    plt.xticks(locs, labels)
    plot.SetupPlot('Temp: '+str(Theta[1,0])+', Bounce: '+str(Theta[3,0]), r'$\theta$ / rad', r'p($\theta)$ / rad$^{-1}$')

def BounceCompare(Temp, Angle, Theta):
    XArr = []
    fl = open('/home/becker/lammps/111/AngDistr.dat', 'r')
    ctr = 0
    for i, line in enumerate(fl):
        a, t, e, nb, x, theta = line.split()
        if i < 35:
            XArr.append(float(x))
        # look at what effect does the bounce threshold have on the ang distr
        if float(t) == Temp and float(a) == Angle:
            Theta[:,ctr] = float(a),float(t),float(e),int(nb),float(x),float(theta)
            ctr += 1


    for num, bnc in enumerate(BncArr):
        PlotArr = []
        for i in range(len(Theta[0,:])):
            if Theta[-1,i] > 0:
                if Theta[3,i] == BncArr[num]:
                    PlotArr.append(Theta[5,i])

        plt.plot(XArr, PlotArr[:size], label=str(bnc)+' Bounces')

    fl.close()
    locs = [0, pi/12, pi/6, pi/4, pi/3, 5*pi/12, pi/2]
    labels = ['0', r'$\frac{\pi}{12}$', r'$\frac{\pi}{6}$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{3}$', r'$\frac{5\pi}{12}$', r'$\frac{\pi}{2}$']
    plt.xticks(locs, labels)
    plot.SetupPlot('Temp: '+str(Theta[1,0])+', Angle: '+str(Theta[0,0]), r'$\theta$ / rad', r'p($\theta)$ / rad$^{-1}$')
