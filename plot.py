# plot.py
import matplotlib.pyplot as plt
import smooth
import math
import stats
from scipy.signal import savgol_filter as sf
import numpy as np
import numpy.ma as ma

def Populations(angle, temp, energy, X, T, Q, C, smflag=1, pltflag=1, nu=5, hlrn=0):
    TIMESTEPS = len(T)
    DT = 1
    if hlrn == 1:
        name = 'a'+angle+'t'+temp+'e'+energy+'HLRN'
    else:
        name = 'a'+angle+'t'+temp+'e'+energy
    if smflag == 1:
        print("Smoothing Populations")
        ctr = 0
        for k in range(0,TIMESTEPS):
            if k % (TIMESTEPS / 10) == 0:
                print(angle, temp, '\t'+str(ctr)+' %')
                ctr += 10
            C[k] = smooth.GaussSmoothing(TIMESTEPS, k, C, DT, nu=nu)
            Q[k] = smooth.GaussSmoothing(TIMESTEPS, k, Q, DT, nu=nu)
            T[k] = smooth.GaussSmoothing(TIMESTEPS, k, T, DT, nu=nu)
        print(angle, temp, '\t100 %')

    if pltflag == 1:
        print('Plot Populations')
        fig, ax = plt.subplots(num='a'+angle+'t'+temp+'e'+energy)
        ax.plot(X*2.5e-2, C, 'r-', label = 'C')
        ax.plot(X*2.5e-2, Q, 'g--', label = 'Q')
        ax.plot(X*2.5e-2, T, 'b-.', label = 'T')
        legend = ax.legend()
        plt.xlabel("time / ps")
        plt.ylabel("")
        plt.title("Angle " + angle + " deg, Energy " + energy + " meV, Temp " + temp + " K")
        plt.tight_layout()
        #plt.savefig('/home/becker/lammps/tmpplot/' + name + 'Population.pdf')
        plt.show(block=True)
        plt.close(fig)

def TransitionPopulations(angle, temp, energy, X, QT, TQ, CT, CQ, TC=[], QC=[], smflag=1, pltflag=1, nu=49, hlrn=0):
    TIMESTEPS = len(QT)
    l = len(QT)
    DT = 1
    fig, ax = plt.subplots(num='a'+angle+'t'+temp+'e'+energy)
    if hlrn == 1:
        name = 'a'+angle+'t'+temp+'e'+energy+'HLRN'
    else:
        name = 'a'+angle+'t'+temp+'e'+energy
    if pltflag == 1:
        ax.plot(X*2.5e-2, QT, 'b--', label = 'QT')
        ax.plot(X*2.5e-2, TQ, 'C1--', label = 'TQ')
        ax.plot(X*2.5e-2, CT, 'g--', label = 'CT')
        ax.plot(X*2.5e-2, CQ, 'r--', label = 'CQ')
    if smflag == 1:
        start = int(5 * l/30)
        Edge = np.arange(-start-5,l-start-5,1)
        print("Smoothing Transition Populations")
        ctr = 0
        # sf = scipy.signal.savgol_filter
        QT[:] = sf(QT, nu, 3, deriv=0)
        TQ[:] = sf(TQ, nu, 3, deriv=0)
        CT[:] = sf(CT, nu, 3, deriv=0)
        CQ[:] = sf(CQ, nu, 3, deriv=0)
        '''
        for k in range(0,TIMESTEPS):# - int(0.01*l)):
            if k % (TIMESTEPS / 10) == 0:
                print(angle, temp, '\t'+str(ctr)+' %')
                ctr += 10

            QT[k] = smooth.GaussSmoothing(TIMESTEPS, k, QT, dt=DT, nu=nu, edge='', X=Edge)
            TQ[k] = smooth.GaussSmoothing(TIMESTEPS, k, TQ, dt=DT, nu=nu, edge='', X=Edge)
            CT[k] = smooth.GaussSmoothing(TIMESTEPS, k, CT, dt=DT, nu=nu, edge='', X=Edge)
            CQ[k] = smooth.GaussSmoothing(TIMESTEPS, k, CQ, dt=DT, nu=nu, edge='', X=Edge)
        print(angle, temp, '\t100 %')
        '''
    # plot the results
    if pltflag == 1:
        print("Plot Transition Populations")
        #fig, ax = plt.subplots(num='a'+angle+'t'+temp+'e'+energy)
        ax.plot(X*2.5e-2, QT, 'b-', label = 'QT')
        ax.plot(X*2.5e-2, TQ, 'C1-', label = 'TQ')
        ax.plot(X*2.5e-2, CT, 'g-', label = 'CT')
        ax.plot(X*2.5e-2, CQ, 'r-', label = 'CQ')
        try:
            ax.plot(X*2.5e-2, TC, '-', label = 'TC')
            ax.plot(X*2.5e-2, QC, '-', label = 'QC')
        except:
            dummy=1

        legend = ax.legend()
        plt.xlabel("time / ps")
        plt.ylabel("")
        plt.title("Angle " + angle + " deg, Energy " + energy + " meV, Temp " + temp + " K")
        #plt.savefig('./check.pdf')
        #if int(temp) == 80:
        #    plt.axis([0,30,0,.1])
        #else:
        #    plt.axis([0,30,0,.3])
        plt.tight_layout()
        #plt.savefig('/home/becker/lammps/tmpplot/' + name + 'Transition.pdf')
        plt.show(block=True)
        plt.close(fig)

    return QT, TQ, CT, CQ, TC, QC


#tpl: # nb, t1, t2, tm, s1, s2, transition, (energy loss ? ) # tuple with information about bounce events
def TransitionRate(angle, temp, energy, X, Ta, Tb, Tc, Td, lblA='', lblB='', lblC='', lblD='', smflag=1, pltflag=1, nu=8, ylbl='', avgflag=0, start=0, end=0, hlrn=0):
    if hlrn == 1:
        name = 'a'+angle+'t'+temp+'e'+energy+'HLRN'
    else:
        name = 'a'+angle+'t'+temp+'e'+energy
    if smflag == 1:
        print("Smoothing Transition Rates")

        for k in range(0,len(Ta)):
            Ta[k] = smooth.GaussSmoothing(len(Ta), k, Ta, dt=1, nu=nu)
            Tb[k] = smooth.GaussSmoothing(len(Ta), k, Tb, dt=1, nu=nu)
            Tc[k] = smooth.GaussSmoothing(len(Ta), k, Tc, dt=1, nu=nu)
            Td[k] = smooth.GaussSmoothing(len(Ta), k, Td, dt=1, nu=nu)

    if pltflag == 1:
        print("Plot Transition Rates")
        fig, ax = plt.subplots(num='a'+angle+'t'+temp+'e'+energy)
        if avgflag == 1:
            avga = stats.AvgRate(start, end, Ta)
            avgb = stats.AvgRate(start, end, Tb)
            avgc = stats.AvgRate(start, end, Tc)
            avgd = stats.AvgRate(start, end, Td)
            Avga = np.array([avga for i in range(len(Ta))])
            Avgb = np.array([avgb for i in range(len(Ta))])
            Avgc = np.array([avgc for i in range(len(Ta))])
            Avgd = np.array([avgd for i in range(len(Ta))])
            ax.plot(X*2.5e-2, Avga, 'b--', label=lblA+'avg')
            ax.plot(X*2.5e-2, Avgb, 'C1--', label=lblB+'avg')
            ax.plot(X*2.5e-2, Avgc, 'g--', label=lblC+'avg')
            ax.plot(X*2.5e-2, Avgd, 'r--', label=lblD+'avg')

        ax.plot(X*2.5e-2, Ta, 'b-', label = lblA)
        ax.plot(X*2.5e-2, Tb, 'C1-', label = lblB)
        ax.plot(X*2.5e-2, Tc, 'g-', label = lblC)
        ax.plot(X*2.5e-2, Td, 'r-', label = lblD)
        legend = ax.legend()
        plt.xlabel("time / ps")
        plt.ylabel(ylbl)
        if hlrn == 1:
            plt.axis([0,60,-0.1,1.5])
        else:
            plt.axis([0,40,-0.1,1.5])
        plt.title("Angle " + angle + " deg, Energy " + energy + " meV, Temp " + temp + " K")
        #plt.savefig('./check.pdf')
        #if int(temp) == 80:
        #    plt.axis([0,30,0,.1])
        #else:
        #    plt.axis([0,30,0,.3])
        plt.tight_layout()
        #plt.savefig('/home/becker/lammps/tmpplot/' + name + 'TransitionRate.pdf')
        plt.show(block=True)
        plt.close(fig)

    return Ta, Tb, Tc, Td

#TODO
def Histogram(emin, emax, nbin, A, B, C, subplts=0, lblA='', lblB='', lblC='', lbl='', Title='', binarr=[], avg=0, std=0,
 totalplt=1, sm_flag=1, pltflag=1, nu=4, edge='none', edgeA='none', edgeB='none', edgeC='none', writeflag=0, flname='', nb=0, action=-1, scl=1):
    # find a better way to create bins, so that 0 meV is always a bin boundary
    bin_len = 3 #TODO
    average = -3000
    stdev = -9000
    if len(binarr) == 0:
        bns = np.linspace(emin, emax, nbin)
    else:
        bns = binarr
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
            sA = np.zeros([len(hTot[0])])
            sB = np.zeros([len(hTot[0])])
            sC = np.zeros([len(hTot[0])])
            for k in range(0,len(hTot[0])):
                sA[k] = smooth.GaussSmoothing(len(hA[0]), k, hA[0], edge=edgeA, nu=nu, X=hTot[1][:-1])
                sB[k] = smooth.GaussSmoothing(len(hB[0]), k, hB[0], edge=edgeB, nu=nu, X=hTot[1][:-1])
                sC[k] = smooth.GaussSmoothing(len(hC[0]), k, hC[0], edge=edgeC, nu=nu, X=hTot[1][:-1])
            if(pltflag == 1):
                if edgeA == 'positive' or edgeA == 'pos':
                        plt.plot(ma.masked_less(hA[1][:-1],0), scl*(sA / norm), label=lblA)
                elif edgeA == 'negative' or edgeA == 'neg':
                        plt.plot(ma.masked_greater(hA[1][:-1],0), scl*(sA / norm), label=lblA)
                else:
                        plt.plot(hA[1][:-1], scl*(sA) / norm, label=lblA)

                if edgeB == 'positive' or edgeB == 'pos':
                        plt.plot(ma.masked_less(hB[1][:-1],0), scl*(sB / norm), label=lblB)
                elif edgeB == 'negative' or edgeB == 'neg':
                        plt.plot(ma.masked_greater(hB[1][:-1],0), scl*(sB / norm), label=lblB)
                else:
                        plt.plot(hB[1][:-1], scl*(sB) / norm, label=lblB)

                if edgeC == 'positive' or edgeC == 'pos':
                        plt.plot(ma.masked_less(hC[1][:-1],0), scl*(sC / norm), label=lblC)
                elif edgeC == 'negative' or edgeC == 'neg':
                        plt.plot(ma.masked_greater(hC[1][:-1],0), scl*(sC / norm), label=lblC)
                else:
                        plt.plot(hC[1][:-1], scl*(sC) / norm, label=lblC)

        else:
            if(pltflag == 1):
                plt.plot(hA[1][:-1], scl*(hA[0] / norm), label=lblA)
                plt.plot(hB[1][:-1], scl*(hB[0] / norm), label=lblB)
                plt.plot(hC[1][:-1], scl*(hC[0] / norm), label=lblC)
    if totalplt == 1:
        if sm_flag == 1:
            sTot = np.zeros([len(hTot[0])])
            for k in range(0, len(hTot[0])):
                sTot[k] = smooth.GaussSmoothing(len(hTot[0]), k, hTot[0], edge=edge, nu=nu, X=hTot[1][:-1])
            if pltflag == 1:
                plt.plot(hTot[1][:-1], scl*sTot / norm, label=lbl)

        else:
            if pltflag == 1:
                plt.plot(hTot[1][:-1], scl*(hTot[0] / norm), label=lbl)
    if pltflag == 1:
        plt.title(Title)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        plt.close()

    if writeflag == 1 and sm_flag == 1 and subplts == 1:
        fl = open(flname,'a')
        fl.write("#action\t bnc\t X\t\t\t total\t\t\t T\t\t\t\t Q\t\t\t\t C\n")
        for i in range(0, len(hTot[0])):
            fl.write("%d\t\t %d\t %e\t %e\t %e\t %e\t %e\n" %(action, nb, hTot[1][i], sTot[i] * scl / norm, sA[i] * scl / norm, sB[i] * scl / norm, sC[i] * scl / norm))

    elif writeflag == 1:
        print('Can only write to file with activated subplots and smoothing')

    if avg == 1:
        average = np.average(hTot[1][:-1], weights=sTot/norm)
    if std == 1 and avg == 1:
        square = np.average(np.square(hTot[1][:-1]), weights=sTot/norm)
        stdev = math.sqrt(square - (average**2)) / nbin
    elif std == 1 and avg != 1:
        print("Error! Need Mean to calculate Std.")

    if subplts == 1:
        if sm_flag ==1:
            return hTot[1], sTot*scl, norm, sA*scl, sB*scl, sC*scl, average, stdev
        else:
            return hTot[1], hTot[0]*scl, norm, hA[0]*scl, hB[0]*scl, hC[0]*scl, average, stdev
    else:
        if sm_flag == 1:
            return hTot[1], sTot*scl, norm, [], [], [], average, stdev
        else:
            return hTot[1], hTot[0]*scl, norm, [], [], [], average, stdev


def EnergyDistribution():
    return 0
