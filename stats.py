import math
import numpy as np

#calculate StdDeviation for Sticking coefficient
def StdDeviation(NQ, NT, N):
    nu = (NQ + NT)# NQ and NT are already normalized to 1; otherwise here would be: (NQ+NT)/N
    return math.sqrt(nu*(1-nu)/N)

def InitSticking(Refl, angle, temp, energy, NUM, printflag=1):
    StickProb = [0,0,0,0]
    ReflCoeff = Refl[-1]
    sigma = StdDeviation(ReflCoeff, 0.0, NUM)
    StickProb[0] = angle
    StickProb[1] = temp
    StickProb[2] = energy
    StickProb[3] = 1.0 - ReflCoeff
    if printflag == 1:
        print(StickProb, '+/-', sigma,'\n')

    return 1.0 - ReflCoeff

def median(Val):
    m = 0
    for i in range(0, len(Val)):
        m += Val[i]
    m /=  len(Val)
    return m

#calculate StdDev for energy development over time
def stdDev(E,avg,NumTraj):
    std = np.zeros([len(E)])
    for e in range(len(E)):
        if E[e] <= 0:
            continue
        std = (E[e] - avg)**2
        std /= NumTraj
    return std
