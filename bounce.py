import numpy as np
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
def Countbounces(Vz, Epar, Enorm, Epot, Bnc):
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
            # nb, tm, s1, t1, E_xy.i, E_z.i, V.i, s2, t2, E_xy.f, E_z.f, V.f, refl_flag (s2 == 1 == C functions as flag)
            Bnc[:,nb] = nb, tm, s1, t1, Epar[t1], Enorm[t1], Epot[t1], s2, t2, Epar[t2], Enorm[t2], Epot[t2], s2
            nb += 1
            i = t2
    Bnc[:,nb] = nb, l+1, s2, l+1, -3000, -3000, -3000, s2, l+1, -3000, -3000, -3000, 1
    return nb

#State: T, Q, C; Trans: QT, CT, TQ, CQ
def Population(nb, bnc1, bnc2, State, Trans, Refl, scale, traj):
    global mistake
    global Tmistake
    global mistakeW
    global flag
    l = len(State[2]) - 1
    assert(bnc2[0] > bnc1[0])
    if int(bnc1[7]) != int(bnc2[2]) and bnc2[1] <= 1199:
        dummyvariable = 1
        #print('check',str(traj), str(nb), str(int(bnc1[7]) - int(bnc2[2])), str(bnc1[1]))
    #assert(bnc2[1] > bnc1[1])
    #if bnc1[7] != bnc2[2]:
        #print("State description not consistent in ", str(traj), str(nb))
    if nb == 0:
        for i in range(0, int(bnc1[1])):
            State[2,i] += 1. * scale
# when particle is in C state, it will not come back to
# any kind of bound state
# therefore we fill the population of C states from the last bounce
# till the end
    if int(bnc1[7]) == 1:# and flag != traj:
        State[2, int(bnc1[1]):len(State[2])] += 1
        #State[2, int(bnc1[1]):int(bnc2[1])] += 1

        if int(bnc1[2]) == -1 and int(bnc1[7]) == 1: #T->C
            Trans[1, int(bnc1[1]):len(State[2])] += 1 * scale

        elif int(bnc1[2]) == 0 and int(bnc1[7]) == 1: #Q->C
            Trans[3, int(bnc1[1]):len(State[2])] += 1. * scale

        elif int(bnc1[2]) == 1 and int(bnc1[7]) == 1: # direct reflection
            Refl[int(bnc1[1]):len(State[2])] += 1. * scale

        #bnc2[:] = 0,0,0,0,0,0,0,0,0,0,0,0,0
        bnc2[:] = nb, l+1, 1, l+1, -3000, -3000, -3000, 1, l+1, -3000, -3000, -3000, 1
        return -3000
        #flag = traj
        #return flag
        #end if
        #return 1
    #end if

    elif int(bnc1[7]) == -1:# and flag != traj:

        if int(bnc1[2]) == 1 and int(bnc1[7]) == -1: #C->T
            if nb > 0:
                mistake += 1
                print("Particle returned to T state.")
                print('traj: ', str(traj), 'nb: ', str(nb), 'err: ', str(mistake))
                #State[2, int(bnc1[1]):len(State[2])] += 1 * scale
                return -6000

            Trans[4, int(bnc1[1]):len(State[0])] += 1 * scale

        elif int(bnc1[2]) == 0 and int(bnc1[7]) == -1: #Q->T
            Trans[2, int(bnc1[1]):len(State[0])] += 1. * scale

        State[0, int(bnc1[1]):int(bnc2[1])] += 1

    elif int(bnc1[7]) == 0:# and flag != traj:

        if int(bnc1[2]) == 1 and int(bnc1[7]) == 0: #C->Q
            if nb > 0:
                mistake += 1
                print("Particle returned to Q state.")
                print('traj: ', str(traj), 'nb: ', str(nb), 'err: ', str(mistake))
                #State[2, int(bnc1[1]):len(State[2])] += 1 * scale
                return -5000

            Trans[5, int(bnc1[1]):len(State[0])] += 1 * scale

        elif int(bnc1[2]) == -1 and int(bnc1[7]) == 0: #T->Q
            Trans[0, int(bnc1[1]):len(State[1])] += 1. * scale

        State[1, int(bnc1[1]):int(bnc2[1])] += 1
    '''
    if int(bnc1[7]) != int(bnc2[2]):
        mistakeW += 1
        print("warning! Diff of States = ", str(bnc1[7] - bnc2[2]))
        print('traj: ', str(traj), 'nb: ', str(nb), 'err: ', str(mistakeW))
        '''
    #if State[a,i]

    #elif int(bnc1[7]) == 1:
    #    print('an error occurred, check bounce.py')
    #    for i in range(int(bnc1[1]), int(bnc2[1])):
    #        State[2,i] += 1 * scale
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

def TransitionRate(dt, Nb, N_ab, T_ab):
    for t in range(0,len(Nb)):
        if Nb[t] == 0:
            T_ab[t] = T_ab[t-1]
        else:
            dq = N_ab[t] / dt
            T_ab[t] = dq / Nb[t]
    return T_ab

# calculate the transition on the grounds of the before and after state
def Transition(s1, s2):
    if s1 == -1 and s2 == 0:
        return 0   # T -> Q: 0
    elif s1 == -1 and s2 == 1:
        return 1   # T -> C: 1
    elif s1 == 0 and s2 == -1:
        return 2   # Q -> T: 2
    elif s1 == 0 and s2 == 1:
        return 3    # Q -> C: 3
    elif s1 == 1 and s2 == -1:
        return 4   # C -> T: 4
    elif s1 == 1 and s2 == 0:
        return 5    # C -> Q: 5
    elif s1 == 1 and s2 == 1:
        return -99
    else:
        return -1

def TransitionValues(num, tpl1, TQ, QT, QC, TC, CT, CQ, scale):
    #deprecated
    if tpl1[3] > len(TQ)-6:
        tpl1[3] = len(TQ)

    if int(tpl1[6]) == 0:
        for i in range(int(tpl1[1]), len(TQ)): # for bounce event, write to array from t1 till t2
            QT[i] += 1 * scale
    if int(tpl1[6]) == 2:
        for i in range(int(tpl1[1]), len(TQ)): # for i in range(int(tpl1[1]), int(tpl1[2])):
            TQ[i] += 1 * scale
    if int(tpl1[6]) == 3:
        for i in range(int(tpl1[1]), len(TQ)):
            CQ[i] += 1 * scale
    if int(tpl1[6]) == 1:
        for i in range(int(tpl1[1]), len(TQ)):
            CT[i] += 1 * scale
    if num == 0:
        return 0
    if int(tpl1[6]) == 4:
        for i in range(int(tpl1[1]), len(TQ)):
            TC[i] += 1 * scale
    if int(tpl1[6]) == 5:
        for i in range(int(tpl1[1]), len(TQ)):
            QC[i] += 1 * scale
    return 1

    '''State = np.zeros([3,TIMESTEPS]) # T, Q, C
    Trans = np.zeros([6,TIMESTEPS]) # QT, CT, TQ, CQ, TC, QC
    TRate = np.zeros([6,TIMESTEPS]) # QT, CT, TQ, CQ, TC, QC
                                       0   1   2   3   4   5
    '''

def IntegratePopulationT(Tnew, State, TRate, dt):
    Tnew[0] = State[0,0]
    for j in range(1, len(State[0])):
        i = j-1
        Tnew[j] = -State[0,i] * (TRate[0][i] + TRate[1][i]) * dt + State[1][i] * TRate[2][i] * dt + State[2][i] * TRate[4][i] * dt# + Tnew[i]
    return Tnew

def IntegratePopulationQ(Qnew, State, TRate, dt):
    Qnew[0] = State[1,0]
    for j in range(1,len(State[1])):
        i = j-1
        Qnew[j] = -State[1][i] * (TRate[3][i] + TRate[2][i]) * dt + State[0][i] * TRate[0][i] * dt + State[2][i] * TRate[5][i] * dt# + Qnew[i]
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
