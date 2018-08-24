import math
import numpy as np
import cfg as g

def Scaling(State, Trans, Refl): #, TRAJ_COUNTER, HLRN_COUNTER, TIMESTEPS_RZ, TIMESTEPS_HLRN):
    for i in range(0, g.TIMESTEPS_HLRN):
        if i < g.TIMESTEPS_RZ:
            for j in range(0,3):
                State[j,i] /= g.TRAJ_COUNTER
            for j in range(0,6):
                Trans[j,i] /= g.TRAJ_COUNTER
            Refl[i] /= g.TRAJ_COUNTER
        else:
            for j in range(0,3):
                State[j,i] /= g.HLRN_COUNTER
            for j in range(0,6):
                Trans[j,i] /= g.HLRN_COUNTER
            Refl[i] /= g.HLRN_COUNTER

#calculate StdDeviation for Sticking coefficient
def StdDeviation(NQ, NT, N):
    nu = (NQ + NT)# NQ and NT are already normalized to 1; otherwise here would be: (NQ+NT)/N
    return math.sqrt(nu*(1-nu)/N)

def InitSticking(InArr1, InArr2, delay_time, angle, temp, energy, NUM, printflag=0, writeflag=0, filename=''):
    StickProb = [0,0,0,0]
    # ReflCoeff = Refl[-1]
    if len(InArr2) == 0:
        ReflCoeff = 1.0 - InArr1[int(delay_time)]
    else:
        ReflCoeff = 1.0 - InArr1[int(delay_time)] + InArr2[int(delay_time)]
    sigma = StdDeviation(ReflCoeff, 0.0, NUM)
    StickProb[0] = angle
    StickProb[1] = temp
    StickProb[2] = energy
    StickProb[3] = 1.0 - ReflCoeff
    if printflag == 1:
        print(StickProb, '+/-', sigma,'\n')
    if writeflag == 1:
        f = open(filename, 'a')
        f.write("%d %d %f %f %f\n" %(angle, temp, energy, 1.0-ReflCoeff, sigma))
        f.close()

    return 1.0 - ReflCoeff, sigma

def median(Val, fltr=0, fltrval=0):
    m = 0
    ctr = 0

    if fltr == 1:
        for i in range(0, len(Val)):
            if Val[i] == fltrval:
                continue
            else:
                m += Val[i]
                ctr += 1
    else:
        for i in range(0, len(Val)):
            m += Val[i]
            ctr += 1
    m /= ctr
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
    l1 = -0.5 * (math.fabs(R[0,0]) + math.fabs(R[1,1]) - math.sqrt( ((math.fabs(R[1,1]) - math.fabs(R[0,0]) ) ** 2)  + 4.* R[0,1] * R[1,0]) )
    l2 = -0.5 * (math.fabs(R[0,0]) + math.fabs(R[1,1]) + math.sqrt( ((math.fabs(R[1,1]) - math.fabs(R[0,0]) ) ** 2)  + 4.* R[0,1] * R[1,0]) )
    return l1, l2

def StdEigenvalue(R, Rstd):
    print("calculate error on eigenvalues")
    a = R[0,0]
    b = R[1,1]
    c = R[0,1]
    d = R[1,0]
    root = math.sqrt( (a-b)**2 + 4*c*d )
    terma = math.fabs(-0.5 - (a-b)/(2. * root)) * Rstd[0,0]
    termb = math.fabs(-0.5 + (a-b)/(2. * root)) * Rstd[1,1]
    termc = math.fabs(-d/root) * Rstd[0,1]
    termd = math.fabs(-c/root) * Rstd[1,0]
    delL = terma + termb + termc + termd
    return delL

def Coeffs(l1, l2, N, R, te):
    c1 = 1. / (l1 -l2) * (N[0][te] - (l2-R[1,1]) * N[1][te] / R[1,0])
    c2 = -1. / (l1 -l2) * (N[0][te] - (l1-R[1,1]) * N[1][te] / R[1,0])
    return c1, c2

def StdCoeffs(l1, l2, N, R, te, lstd, Rstd, total):
    print("calculate error on coefficients")
    denom = (l1-l2)
    delN0 = StdDeviation(N[0][te], 0, total)
    delN1 = StdDeviation(N[1][te], 0, total)
    term_l1 = math.fabs(lstd / denom**2 * (N[0][te] - (l2-R[1,1]) * N[1][te]/R[1,0]))
    term_l2 = math.fabs(lstd * (N[1][te]*l1 - N[1][te]*R[1,1] - N[0][te]*R[1,0]) / (R[1,0] * denom**2))
    term_N0 = math.fabs(delN0 / denom)
    term_N1 = math.fabs(delN1 * (R[1,1] - l2) / (R[1,0] * denom))
    term_R11 = math.fabs(Rstd[1,1] * N[1][te] / (R[1,0] * denom))
    term_R10 = math.fabs(Rstd[1,0] * N[1][te] * (l2 - R[1,1]) / (R[1,0]**2 * denom))
    delC1 = term_l1 + term_l2 + term_N0 + term_N1 + term_R11 + term_R10

    term_l1 = math.fabs(lstd * (N[1][te]*R[1,1] + N[0][te]*R[1,0] - l2*N[1][te]) / (R[1,0] * denom))
    term_l2 = math.fabs(lstd * (N[0][te] - N[1][te]*(l1-R[1,1])/R[1,0]) / denom**2)
    term_N0 = math.fabs(delN0 / denom)
    term_N1 = math.fabs(delN1 * (l1 - R[1,1]) / (R[1,0] * denom))
    term_R11 = math.fabs(Rstd[1,1] * N[1][te] / (R[1,0] * denom))
    term_R10 = math.fabs(Rstd[1,0] * N[1][te] * (l1 - R[1,1]) / (R[1,0]**2 * denom))
    delC2 = term_l1 + term_l2 + term_N0 + term_N1 + term_R11 + term_R10
    return delC1, delC2

def PopN(c1, c2, l1, l2, R, N, Max, te):
    '''for t in range(0,Max):
        # check for conversion of t to right units!
        N[0][t] = c1 * (l1 - R[1,1]) * math.exp(l1 * t) + c2 * (l2 - R[1,1]) * math.exp(l2 * t)
        '''
    # alternatively:
    tarr = np.arange(0,Max) - te
    tarr = tarr * 0.025  / 6.53
    N[0] = c1 * (l1 - R[1,1]) * np.exp(l1 * tarr) + c2 * (l2 - R[1,1]) * np.exp(l2 * tarr)
    N[1] = c1 * R[1,0] * np.exp(l1 * tarr) + c2 * R[1,0] * np.exp(l2 * tarr)
    return N


def AvgRate(start, end, TRate):
    avg = 0
    M = end - start
    for i in range(start, end):
        avg += TRate[i] / M
    return avg

def StdRate(start, end, TRate, avg):
    std = 0
    M = end - start

    for i in range(start, end):
        std += (TRate[i] - avg) ** 2 / (M-1)
    return np.sqrt(std) / np.sqrt(M)                # effectively we calculate sigma/sqrt(M)
                                                    # standard error of the mean

def CalcR(start, end, R, Rstd, TRate):
    avgt0 = AvgRate(start, end, TRate[0])
    avgt1 = AvgRate(start, end, TRate[1])
    R[0,0] = -(avgt1 + avgt0)
    Rstd[0,0] = StdRate(start, end, TRate[1], avgt1) + StdRate(start, end, TRate[0], avgt0)

    avgt2 = AvgRate(start, end, TRate[2])
    avgt3 = AvgRate(start, end, TRate[3])
    R[1,1] = -(avgt3 + avgt2)
    Rstd[1,1] = StdRate(start, end, TRate[3], avgt3) + StdRate(start, end, TRate[2], avgt2)

    R[0,1] = avgt2
    Rstd[0,1] = StdRate(start, end, TRate[2], avgt2)
    R[1,0] = avgt0
    Rstd[1,0] = StdRate(start, end, TRate[0], avgt0)
    return R

def ResTime(l1, N00, N10, R):
    Delta = R[0,1] * R[1,0] / math.fabs(R[1,1])
    deltalambda = math.fabs(R[1,1]) - math.fabs(R[0,0]) + 2. * Delta
    gamma = (1./math.e)
    N0 = (1. - gamma) * ( N00 + R[0,1] / math.fabs(R[1,1]) * N10 ) + (N10 + R[1,0] * N00 / deltalambda)
    tR = 1./l1 * math.log(gamma * (N00 + N10) / N0)
    return tR * 6.53

def StdResTime(l1, delL, N00, N10, R):
    #TODO: error on N0
    Delta = R[0,1] * R[1,0] / math.fabs(R[1,1])
    gamma = (1./math.e)
    deltalambda = math.fabs(R[1,1]) - math.fabs(R[0,0]) + 2. * Delta
    N0 = (1. - gamma) * ( N00 + R[0,1] / math.fabs(R[1,1]) * N10 ) + (N10 + R[1,0] * N00 / deltalambda)
    std = delL / l1**2 * math.log(gamma * (N00 + N10) / N0)
    return math.fabs(std * 6.53)

def CalcSolution(fp, State, start, end, R, Rstd, TRate, num, multiplier=2):
    N = np.zeros([2, int(g.TIMESTEPS*multiplier)])
    N[0][fp] = State[0][fp]
    N[1][fp] = State[1][fp]
    Time2 = np.arange(0, int(g.TIMESTEPS*multiplier))
    CalcR(start, end, R[:,:,0], Rstd[:,:,0], TRate[:, :num])
    l1, l2 = Eigenvalue(R[:,:,0])
    delL = StdEigenvalue(R[:,:,0], Rstd[:,:,0])
    c1, c2 = Coeffs(l1, l2, N, R[:,:,0], fp)
    dc1, dc2 = StdCoeffs(l1, l2, N, R[:,:,0], fp, delL, Rstd[:,:,0], g.TRAJ_COUNTER)
    N = PopN(c1, c2, l1, l2, R[:,:,0], N, int(g.TIMESTEPS*multiplier), fp)

    return N, l1, l2, delL, c1, dc1, c2, dc2, Time2

def TMAC(vxi, vyi, vxf, vyf):
    tmac = 1. - (vxi * vxf + vyi * vyf) / (vxi**2 + vyi**2)
    return tmac

def nrgCorr(X, Y):
    cov = E(X*Y) - E(X)*E(Y)
    varX = var(X)
    varY = var(Y)
    corr = cov / (np.sqrt(varX * varY) )
    return corr

def TMACstd(tmac, avgtmac):
    ctr = 0
    sigma = 0
    for i in range(len(tmac)):
        if tmac[i] == 0:
            continue
        sigma += (tmac[i] - avgtmac)**2
        ctr += 1
    sigma = np.sqrt(1/(ctr-1) * sigma)
    stdtmac = sigma / np.sqrt(ctr)
    return stdtmac
