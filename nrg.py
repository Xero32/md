import numpy as np
import smooth
import matplotlib.pyplot as plt

def EnergyLoss(nb, delT, delQ, delC, length, Ei, Ef, state2, flag):  # implementation: E == NRG[3:8,:,:]
    t = q = c = 0
    for j in range(0,length):
        if int(state2[nb,j]) == -1:# and int(flag[nb,j]) == 0:
            delT[t] = (Ef[0,nb,j]+Ef[1,nb,j]+Ef[2,nb,j])*1e3 - (Ei[0,nb,j]+Ei[1,nb,j]+Ei[2,nb,j])*1e3 # 1e3 conversion from eV to meV
            if delT[t] == 0:
                #print("Delta T energy = 0", str(t), str(j))
                continue
            t += 1
        elif int(state2[nb,j]) == 0:# and int(flag[nb,j]) == 0:
            delQ[q] = (Ef[0,nb,j]+Ef[1,nb,j]+Ef[2,nb,j])*1e3 - (Ei[0,nb,j]+Ei[1,nb,j]+Ei[2,nb,j])*1e3 # 1e3 conversion from eV to meV
            if delQ[q] == 0:
                #print("Delta Q energy = 0", str(q), str(j))
                continue
            q += 1
        elif int(state2[nb,j]) == 1:# and int(flag[nb,j]) == 0:
            delC[c] = (Ef[0,nb,j]+Ef[1,nb,j]+Ef[2,nb,j])*1e3 - (Ei[0,nb,j]+Ei[1,nb,j]+Ei[2,nb,j])*1e3 # 1e3 conversion from eV to meV
            if delC[c] == 0:
                continue
            c += 1


    return t, q, c

def paraEnergyDistr(nb, E_T, E_Q, E_C, length, E, state2, flag):
    t = q = c = 0
    for j in range(0,length):
        if (int(state2[nb,j])) == -1:# and int(flag[nb,j]) == 0:
            E_T[t] = E[0,nb,j]*1e3
            t += 1
        elif (int(state2[nb,j])) == 0:# and int(flag[nb,j]) == 0:
            E_Q[q] = E[0,nb,j]*1e3
            if E_Q[q] == 0:
                #print("Parallel Q energy = 0", str(q), str(j))
                continue
            q += 1
        elif (int(state2[nb,j])) == 1:# and int(flag[nb,j]) == 0:
            E_C[c] = E[0,nb,j]*1e3
            c += 1
    return t, q, c

def normEnergyDistr(nb, E_T, E_Q, E_C, length, E, state2, flag):
    t = q = c = 0
    for j in range(0,length):
        if (int(state2[nb,j])) == -1:# and int(flag[nb,j]) == 0:
            E_T[t] = E[1,nb,j]*1e3
            t += 1
        elif (int(state2[nb,j])) == 0:# and int(flag[nb,j]) == 0:
            E_Q[q] = E[1,nb,j]*1e3
            if E_Q[q] == 0:
                #print("Normal Q energy = 0", str(q), str(j))
                continue
            q += 1
        elif (int(state2[nb,j])) == 1:# and int(flag[nb,j]) == 0:
            E_C[c] = E[1,nb,j]*1e3
            if E_C[c] == 0:
                #print("Normal C energy = 0", str(q), str(j))
                continue
            c += 1
    return t, q, c

def potEnergyDistr(nb, E_T, E_Q, E_C, length, E, state2, flag):
    t = q = c = 0
    for j in range(0,length):
        if (int(state2[nb,j])) == -1:# and int(flag[nb,j]) == 0:
            E_T[t] = E[2,nb,j]*1e3
            t += 1
        elif (int(state2[nb,j])) == 0:# and int(flag[nb,j]) == 0:
            E_Q[q] = E[2,nb,j]*1e3
            if E_Q[q] == 0:
                #print("Potential Q energy = 0", str(q), str(j))
                continue
            q += 1
        if (int(state2[nb,j])) == 1:# and int(flag[nb,j]) == 0:
            E_C[c] = E[2,nb,j]*1e3
            c += 1

    return t, q, c

def totEnergyDistr(nb, E_T, E_Q, E_C, length, E, state2, flag):
    t = q = c = 0

    for j in range(0,length):
        if (int(state2[nb,j])) == -1:# and int(flag[nb,j]) == 0:
            E_T[t] = (E[0,nb,j]+E[1,nb,j]+E[2,nb,j])*1e3
            t += 1
        elif (int(state2[nb,j])) == 0:# and int(flag[nb,j]) == 0:
            E_Q[q] = (E[0,nb,j]+E[1,nb,j]+E[2,nb,j])*1e3
            if E_Q[q] == 0:
                #print("Total Q energy = 0", str(q), str(j))
                continue
            q += 1
        if (int(state2[nb,j])) == 1:# and int(flag[nb,j]) == 0:
            E_C[c] = (E[0,nb,j]+E[1,nb,j]+E[2,nb,j])*1e3
            c += 1

    return t, q, c

def PlotHist(emin, emax, nbin, A, B, C, subplts=0, lblA='', lblB='', lblC='', lbl='', Title='', totalplt=1, sm_flag=1, sm_win='hanning', win_len=7, pltflag=1):
    binwidth = 220.0/40.0
    win_end = win_len - 3
    # find a better way to create bins, so that 0 meV is always a bin boundary
    bns = np.linspace(emin, emax, nbin)
    Tot = np.concatenate((A,B,C))
    hTot = np.histogram(Tot, bins=bns)
    hTotNorm = np.histogram(Tot, density=True, bins=bns)
    #norm = hTot[0][nbin//2] / hTotNorm[0][nbin//2]
    norm = max(hTot[0]) / max(hTotNorm[0])
    if subplts == 1:
        hA = np.histogram(A, bins=bns)
        hB = np.histogram(B, bins=bns)
        hC = np.histogram(C, bins=bns)
        if sm_flag == 1:
            sA = smooth.smoothing(hA[0] / norm, window_len=win_len, window=sm_win)
            sB = smooth.smoothing(hB[0] / norm, window_len=win_len, window=sm_win)
            sC = smooth.smoothing(hC[0] / norm, window_len=win_len, window=sm_win)
            if(pltflag == 1):
                plt.plot(hA[1][:-1], 1e1*(sA[2:-win_end]), label=lblA)
                plt.plot(hB[1][:-1], 1e1*(sB[2:-win_end]), label=lblB)
                plt.plot(hC[1][:-1], 1e1*(sC[2:-win_end]) ,label=lblC)
        else:
            if(pltflag == 1):
                plt.plot(hA[1][:-1], 1e1*(hA[0] / norm), label=lblA)
                plt.plot(hB[1][:-1], 1e1*(hB[0] / norm), label=lblB)
                plt.plot(hC[1][:-1], 1e1*(hC[0] / norm),label=lblC)
    if totalplt == 1:
        if sm_flag == 1:
            sTot = smooth.smoothing(hTot[0] / norm, window_len=win_len, window=sm_win)
            if pltflag == 1:
                plt.plot(hTot[1][:-1], 1e1*(sTot[2:-win_end]), label=lbl)
        else:
            if pltflag == 1:
                plt.plot(hTot[1][:-1], 1e1*(hTot[0] / norm), label=lbl)
    if pltflag == 1:
        plt.title(Title)
        plt.legend()
        plt.grid(True)
        plt.show()
        plt.close()

    if subplts == 1:
        if sm_flag ==1:
            return hTot[1], sTot[2:-win_end]*1e1, norm, sA[2:-win_end]*1e1, sB[2:-win_end]*1e1, sC[2:-win_end]*1e1
        else:
            return hTot[1], hTot[0]*1e1, norm, hA[0]*1e1, hB[0]*1e1, hC[0]*1e1
    else:
        if sm_flag == 1:
            return hTot[1], sTot[2:-win_end]*1e1, norm, [], [], []
        else:
            return hTot[1], hTot[0]*1e1, norm, [], [], []
"""
def PlotAnim(Efunc, hA=[], hB=[], hC=[], hTot, ):
    fig, ax = plt.subplots()


    Efunc(nb, A, B, C, length, E, state2, flag)
"""
