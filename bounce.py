import numpy as np
import math
from stats import median

# TODO
# calculate the transition-histogram for Plotting
#

flag = 0
e0 = 1.60217662e-19 #(in Coulomb)
amu = 1.66054e-27 #(amu to kg)
e0inv = 1 / e0
mistake = 0
Tmistake = 0
mistakeW = 0
flag = -1
counter = 0

def normalKE(v):
    return 0.5 * 40 * amu * v * 100 * v * 100 * e0inv
def parallelKE(vvx, vvy):
    return 0.5 * 40 * amu * e0inv * (vvx*100*vvx*100 + vvy*100*vvy*100)
def normalEnergySum(nke, pe):
    return nke + pe
def totalKE(nke, pke):
    return nke + pke



def State(nke,pke,pot):
    if (nke+pot) > 0:
        return 1
    else:
        if (nke+pke+pot) >= 0:
            return 0
        else:
            return -1

# nb, tm, s1, t1, E_xy.i, E_z.i, V.i, s2, t2, E_xy.f, E_z.f, V.f, refl_flag
def Countbounces(Vz, Epar, Enorm, Epot, Bnc, current):
    global counter
    l = Vz.size - 1
    nb = 0
    t1 = t2 = tm = 0
    vmin = vmax = 0
    for i in range(0,l):
        if i <= t2:
            continue
        if(Vz[i] < 0 and Vz[i+1] > 0):
            vmin = Vz[i]
            vmax = Vz[i+1]

            #make sure we stay in the array
            if i < 20:
                lbound = 0
            else:
                lbound = i-20

            if i >= l-20:
                rbound = l
            else:
                rbound = i + 20

            #determine min and max velocities
            for j in range(i, lbound, -1):
                if Vz[j] < vmin:
                    vmin = Vz[j]
                    t1 = j
            for j in range(i, rbound):
                if Vz[j] > vmax:
                    vmax = Vz[j]
                    t2 = j

            tm = int(0.5*(t2+t1))
            s1 = State(Enorm[t1], Epar[t1], Epot[t1])
            s2 = State(Enorm[t2], Epar[t2], Epot[t2])
            #print(str(t1), str(t2))
            if nb > 0:
                if s1 == 1 and s2 != 1:
                    counter += 1
                    print('Miscounted state. Getting rid of trajectory information. ' + str(current) + ' ' + str(counter))
                    Bnc[:,nb] = -9999, -100, -100, -100, -3000, -3000, -3000, -100, -100, -3000, -3000, -3000, -9999
                    return -9999
            # nb, tm, s1, t1, E_xy.i, E_z.i, V.i, s2, t2, E_xy.f, E_z.f, V.f, refl_flag (s2 == 1 == C functions as flag)
            Bnc[:,nb] = nb, tm, s1, t1, Epar[t1], Enorm[t1], Epot[t1], s2, t2, Epar[t2], Enorm[t2], Epot[t2], s2
            nb += 1
            i = t2
    Bnc[:,nb] = nb, l+1, s2, l+1, -3000, -3000, -3000, s2, l+1, -3000, -3000, -3000, 1
    return nb

#State: T, Q, C; Trans: QT, CT, TQ, CQ
def Population(nb, bnc0, bnc1, bnc2, State, Trans, Refl, scale, traj):
    global mistake
    global Tmistake
    global mistakeW
    global flag
    l = len(State[2]) - 1
    assert(bnc2[0] > bnc1[0])


    if bnc1[7] != bnc2[2]:
        dummy = 1
        #print(str(traj), str(nb))
        #print('Error!')

    if nb == 0:
        for i in range(0, int(bnc1[1])):
            State[2,i] += 1. * scale

        if int(bnc1[7]) == 1:
            State[2, int(bnc1[1]):len(State[2])] += 1 * scale
            Refl[int(bnc1[1]):len(State[2])] += 1 * scale
            return -3000
        elif int(bnc1[7]) == 0:
            State[1, int(bnc1[1]):int(bnc2[1])] += 1 * scale
            Trans[5, int(bnc1[1]):len(State[1])] += 1 * scale
        elif int(bnc1[7]) == -1:
            State[0, int(bnc1[1]):int(bnc2[1])] += 1 * scale
            Trans[4, int(bnc1[1]):len(State[2])] += 1 * scale
# when particle is in C state, it will not come back to
# any kind of bound state
# therefore we fill the population of C states from the last bounce
# till the end
    else:
        if int(bnc1[7]) == 1:# and flag != traj:
            State[2, int(bnc1[1]):len(State[2])] += 1 * scale
            #State[2, int(bnc1[1]):int(bnc2[1])] += 1

            #TODO
            #if int(bnc1[2]) == -1 and int(bnc1[7]) == 1: #T->C
            if int(bnc0[7]) == -1 and int(bnc1[7]) == 1: #T->C
                Trans[1, int(bnc1[1]):len(State[2])] += 1 * scale
            #TODO
            #elif int(bnc1[2]) == 0 and int(bnc1[7]) == 1: #Q->C
            elif int(bnc0[7]) == 0 and int(bnc1[7]) == 1: #Q->C
                Trans[3, int(bnc1[1]):len(State[2])] += 1. * scale
                #print('it happened')
            #TODO
            #elif int(bnc1[2]) == 1 and int(bnc1[7]) == 1: # direct reflection
            elif int(bnc0[7]) == 1 and int(bnc1[7]) == 1: # direct reflection
                print('supposed direct reflection at later bounce', str(nb), str(traj))
                Refl[int(bnc1[1]):len(State[2])] += 1. * scale

            bnc2[:] = nb, l+1, 1, l+1, -3000, -3000, -3000, 1, l+1, -3000, -3000, -3000, 1
            return -3000
            #end if
        #end if

        elif int(bnc1[7]) == -1:# and flag != traj:

            if int(bnc1[2]) == 1 and int(bnc1[7]) == -1: #C->T
                if nb > 0:
                    mistake += 1
                    print("Particle returned to T state.")
                    print('traj: ', str(traj), 'nb: ', str(nb), 'err: ', str(mistake))
                    #State[2, int(bnc1[1]):len(State[2])] += 1 * scale
                    Trans[4, int(bnc1[1])+1:len(State[0])] -= 1 * scale
                    return -6000

                Trans[4, int(bnc1[1]):len(State[0])] += 1 * scale
            #TODO
            #elif int(bnc1[2]) == 0 and int(bnc1[7]) == -1: #Q->T
            elif int(bnc0[7]) == 0 and int(bnc1[7]) == -1: #Q->T
                Trans[2, int(bnc1[1]):len(State[0])] += 1. * scale

            State[0, int(bnc1[1]):int(bnc2[1])] += 1 * scale

        elif int(bnc1[7]) == 0:# and flag != traj:

            #TODO
            #if int(bnc1[2]) == 1 and int(bnc1[7]) == 0: #C->Q
            if int(bnc0[7]) == 1 and int(bnc1[7]) == 0: #C->Q
                if nb > 0:
                    mistake += 1
                    print("Particle returned to Q state.")
                    print('traj: ', str(traj), 'nb: ', str(nb), 'err: ', str(mistake))
                    #State[2, int(bnc1[1]):len(State[2])] += 1 * scale
                    Trans[5, int(bnc1[1])+1:len(State[0])] -= 1 * scale
                    return -5000

                Trans[5, int(bnc1[1]):len(State[0])] += 1 * scale
            #TODO
            #elif int(bnc1[2]) == -1 and int(bnc1[7]) == 0: #T->Q
            elif int(bnc0[7]) == -1 and int(bnc1[7]) == 0: #T->Q
                Trans[0, int(bnc1[1]):len(State[1])] += 1. * scale

            State[1, int(bnc1[1]):int(bnc2[1])] += 1 * scale

    return flag


# nb, tm, s1, t1, E_xy.i, E_z.i, V.i, s2, t2, E_xy.f, E_z.f, V.f, refl_flag
def Evalbounce(nb, Bnc, Epar, Enorm, Epot):
    t1 = int(Bnc[3,nb])
    t2 = int(Bnc[8,nb])
    tm = int(Bnc[1,nb])
    s1 = State(Enorm[t1], Epar[t1], Epot[t1])
    s2 = State(Enorm[t2], Epar[t2], Epot[t2])
    Bnc[2,nb] = s1
    Bnc[7,nb] = s2
    return s2

# try and find the diff quotient for dt, which reaches from one bounce to the next one
#Deprecated
def TransitionRate(dt, Nb, N_ab, T_ab):
    for t in range(0,len(Nb),dt):
        if Nb[t] == 0:
            T_ab[t] = T_ab[t-1]
        else:
            dq = N_ab[t] / dt
            T_ab[t] = dq / Nb[t]
    return T_ab


# TODO
# numerical differentiation, where we iteratively look for the most recent change in slope
def DifferentiateT(dt, Nb, N_ab, T_ab):
    ddt = dt * 0.025    # [ps]
    ddt = ddt / 6.53    # [t_0]
    # this constitutes unit conversion to t_0 = 6.53 ps
    # effectively T * tt * gamma = T * t_0
    # where gamma = 6.53/0.025 and tt = 0.025

    h = dt
    ctr = 2
    start = 10 //0.025
    a = 0
    b = 0
    save = 0
    for t in range(int(start), len(Nb)):
        if Nb[t] == 0:
            continue
        if N_ab[t] == 0:
            continue
        else:
            try:
                a = t+dt
                #diff = N_ab[a] - N_ab[save]
                b = t-dt
                diff = N_ab[a] - N_ab[b]

                while diff == 0:
                    right = N_ab[b-dt]
                    diff = N_ab[t+dt] - right
                    b -= dt
                    save = b
                    ctr += 1
                    #print(str(diff), str(t), str(b))

                T_ab[t] = diff / (ctr * ddt * Nb[t])
                ctr = 2
            except:
                a = len(N_ab) - 1
                b = t-dt
                diff = N_ab[a] - N_ab[b]

                while diff == 0:
                    right = N_ab[b-dt]
                    diff = N_ab[a] - right
                    b -= dt
                    save = b
                    ctr += 1
                    #print(print(str(diff)), str(t), str(b), str(a))

                T_ab[t] = diff / (ctr * ddt * Nb[t])
                ctr = 2
    return T_ab

def AlternTransition(Trans, State, fixpoint, TRateA, Max):
    for t in range(0,Max):
        diff = Trans[t] - Trans[fixpoint]
        dt = t - fixpoint
        ddt = dt * 0.025 / 6.53
        if ddt != 0 and State[t] != 0:
            TRateA[t] = diff / State[t] / ddt
        else:
            TRateA[t] = 0
    return TRateA

def TransitionRateE(Trans, State, te, TRateE, Max):
    for t in range(0, Max):
        if((State[t] - State[int(te)]) == 0):
            continue
        TRateE[t] = (Trans[t] - Trans[int(te)])/ math.fabs(State[t] - State[int(te)]) / 6.53
    return TRateE

    '''State = np.zeros([3,TIMESTEPS]) # T, Q, C
    Trans = np.zeros([6,TIMESTEPS]) # QT, CT, TQ, CQ, TC, QC
    TRate = np.zeros([6,TIMESTEPS]) # QT, CT, TQ, CQ, TC, QC
                                       0   1   2   3   4   5
    '''

def IntegratePopulationT(Tnew, State, TRate, dt):
    Tnew[0] = State[0,0]
    for j in range(1, len(State[0])):
        i = j-1
        Tnew[j] = -State[0,i] * (TRate[0][i] + TRate[1][i]) * dt + State[1][i] * TRate[2][i] * dt# + State[2][i] * TRate[4][i] * dt# + Tnew[i]
    return Tnew

def IntegratePopulationQ(Qnew, State, TRate, dt):
    Qnew[0] = State[1,0]
    for j in range(1,len(State[1])):
        i = j-1
        Qnew[j] = -State[1][i] * (TRate[3][i] + TRate[2][i]) * dt + State[0][i] * TRate[0][i] * dt# + State[2][i] * TRate[5][i] * dt# + Qnew[i]
    return Qnew

def IntegratePopulationC(Cnew, State, TRate, dt):
    Cnew[0] = State[2,0]
    for j in range(1, len(State[2])):
        i = j-1
        Cnew[j] = 1 + TRate[1][i] * State[0][i] * dt + TRate[3][i] * State[1][i] * dt - State[2][i] * (TRate[4][i] + TRate[5][i]) * dt# + Cnew[i]
    return Cnew

def IntegrationStep(told, T, Q, T_QT, T_CT, T_TQ, dt):
    tnew = -T * (T_QT + T_CT) * dt + Q * T_TQ * dt + told
    return tnew

def IntegrateTransition(ABnew, T_AB, B, dt):
    for j in range(1, len(B)):
        i = j-1
        ABnew[j] = T_AB[i] * B[i] * dt
    return ABnew

def Integrate(Anew, A, dt):
    for j in range(1, len(A)):
        i=j-1
        Anew[j] = A[i] * dt + Anew[i]
    return Anew

#calculate average bounce time + stddev
def BncTime(Bnc, ntraj, nbnc, steps, Arr):
    avg = 0
    ctr = 0
    for n in range(0,nbnc):
        for i in range(0,ntraj):
            if(Bnc[1,n,i] > 0 and Bnc[1,n,i] < steps):
                avg += Bnc[1,n,i]
                ctr += 1
            #end if(count real bounces)
        #end for(scan each trajectory)
        Arr[0,n] = avg / ctr
        Arr[2,n] = ctr
        print(ctr)
        ctr = 0
        avg = 0
    #end for(scan each bounce)

    for n in range(0,nbnc):
        for i in range(0,ntraj):
            if(Bnc[1,n,i] > 0 and Bnc[1,n,i] < steps):
                Arr[1,n] += (Bnc[1,n,i] - Arr[0,n])**2
        Arr[1,n] = np.sqrt(Arr[1,n] / (Arr[2,n]-1))

    return Arr
