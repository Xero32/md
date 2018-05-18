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
        std = np.sqrt((E[e] - avg)**2 / NumTraj)
    return std

def Eigenvalue(R):
    l1 = -0.5 * (math.fabs(R[0,0]) + math.fabs(R[1,1]) - math.sqrt( (math.fabs(R[1,1]) - math.fabs(R[0,0]) ** 2 ) + 4.* R[0,1] * R[1,0]) )
    l2 = -0.5 * (math.fabs(R[0,0]) + math.fabs(R[1,1]) + math.sqrt( (math.fabs(R[1,1]) - math.fabs(R[0,0]) ** 2 ) + 4.* R[0,1] * R[1,0]) )
    return l1, l2

def Coeffs(l1, l2, N, R, te):
    c1 = 1. / (l1 -l2) * (N[0][te] - (l2-R[1,1]) * N[1][te] / R[1,0])
    c2 = -1. / (l1 -l2) * (N[0][te] - (l2-R[1,1]) * N[1][te] / R[1,0])
    return c1, c2

def PopN(c1, c2, l1, l2, R, Max):
    '''for t in range(0,Max):
        # check for conversion of t to right units!
        N[0][t] = c1 * (l1 - R[1,1]) * math.exp(l1 * t) + c2 * (l2 - R[1,1]) * math.exp(l2 * t)
        '''
    # alternatively:
    tarr = np.arange(0,Max) * 0.025 # / 6.53
    N[0] = c1 * (l1 - R[1,1]) * np.exp(l1 * tarr) + c2 * (l2 - R[1,1]) * np.exp(l2 * tarr)
    N[1] = c1 * R[1,0] * np.exp(l1 * tarr) + c2 * R[1,0] * exp(l2 * tarr)


def AvgRate(start, end, TRate):
    avg = 0
    M = end - start
    for i in range(start, end):
        avg += TRate[i] / M

def StdRate(start, end, TRate, avg):
    std = 0
    M = end - start
    for i in range(start, end):
        std = (TRate[i] - avg) ** 2
        std /= M - 1

def CalcR(start, end, R, Rstd, TRate):
    R[0,0] = -stats.AvgRate(start, end, TRate[1]) - stats.AvgRate(start, end, Trate[0])
    Rstd[0,0] = stats.StdRate(start, end, TRate[1], R[0,0]) + stats.StdRate(start, end, TRate[0], R[0,0])
    R[1,1] = -stats.AvgRate(start, end, TRate[3]) - stats.AvgRate(start, end, Trate[2])
    Rstd[1,1] = stats.StdRate(start, end, TRate[3], R[1,1]) + stats.StdRate(start, end, TRate[2], R[1,1])
    R[0,1] = -stats.AvgRate(start, end, TRate[2])
    Rstd[0,1] = stats.StdRate(start, end, TRate[2], R[0,1])
    R[1,0] = stats.AvgRate(start, end, TRate[0])
    Rstd[1,0] = stats.stdRate(start, end, TRate[0], R[1,0])
    return R
