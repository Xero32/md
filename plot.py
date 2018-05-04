# plot.py
import matplotlib.pyplot as plt
import smooth
import numpy as np

def Populations(angle, temp, energy, X, T, Q, C, smflag=1, pltflag=1, nu=2):
    TIMESTEPS = 1200
    DT = 1
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
        #plt.savefig('/home/becker/lammps/111/' + name + 'BounceSmooth.pdf')
        plt.show(block=True)
        plt.close(fig)

def TransitionPopulations(angle, temp, energy, X, QT, TQ, CT, CQ, TC, QC, smflag=1, pltflag=1, nu=2):
    TIMESTEPS = 1200
    DT = 1
    if smflag == 1:
        print("Smoothing Transition Populations")
        ctr = 0
        for k in range(0,TIMESTEPS):
            if k % (TIMESTEPS / 10) == 0:
                print(angle, temp, '\t'+str(ctr)+' %')
                ctr += 10
            QT[k] = smooth.GaussSmoothing(TIMESTEPS, k, QT, DT, nu=nu)
            TQ[k] = smooth.GaussSmoothing(TIMESTEPS, k, TQ, DT, nu=nu)
            CT[k] = smooth.GaussSmoothing(TIMESTEPS, k, CT, DT, nu=nu)
            CQ[k] = smooth.GaussSmoothing(TIMESTEPS, k, CQ, DT, nu=nu)
        print(angle, temp, '\t100 %')
    # plot the results
    if pltflag == 1:
        print("Plot Transition Populations")
        fig, ax = plt.subplots(num='a'+angle+'t'+temp+'e'+energy)
        ax.plot(X*2.5e-2, QT, '-', label = 'QT')
        ax.plot(X*2.5e-2, TQ, '-', label = 'TQ')
        ax.plot(X*2.5e-2, CT, '-', label = 'CT')
        ax.plot(X*2.5e-2, CQ, '-', label = 'CQ')
        ax.plot(X*2.5e-2, TC, '-', label = 'TC')
        ax.plot(X*2.5e-2, QC, '-', label = 'QC')
        legend = ax.legend()
        plt.xlabel("time / ps")
        plt.ylabel("")
        plt.title("Angle " + angle + " deg, Energy " + energy + " meV, Temp " + temp + " K")
        #plt.savefig('./check.pdf')
        #if int(temp) == 80:
        #    plt.axis([0,30,0,.1])
        #else:
        #    plt.axis([0,30,0,.3])
        #plt.savefig('/home/becker/lammps/111/' + name + 'Transition.pdf')
        plt.show(block=True)
        plt.close(fig)


#tpl: # nb, t1, t2, tm, s1, s2, transition, (energy loss ? ) # tuple with information about bounce events
def TransitionRate(angle, temp, energy, X, Ta, Tb, Tc, Td, Te, Tf, lblA='', lblB='', lblC='', lblD='', smflag=1, pltflag=1, nu=10):
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
        ax.plot(X*2.5e-2, Ta, '-', label = lblA)
        ax.plot(X*2.5e-2, Tb, '-', label = lblB)
        ax.plot(X*2.5e-2, Tc, '-', label = lblC)
        ax.plot(X*2.5e-2, Td, '-', label = lblD)
        ax.plot(X*2.5e-2, Te, '-', label = 'T_TC')
        ax.plot(X*2.5e-2, Tf, '-', label = 'T_QC')
        legend = ax.legend()
        plt.xlabel("time / ps")
        plt.ylabel("")
        plt.title("Angle " + angle + " deg, Energy " + energy + " meV, Temp " + temp + " K")
        #plt.savefig('./check.pdf')
        #if int(temp) == 80:
        #    plt.axis([0,30,0,.1])
        #else:
        #    plt.axis([0,30,0,.3])
        #plt.savefig('/home/becker/lammps/111/' + name + 'TransitionRate.pdf')
        plt.show(block=True)
        plt.close(fig)

#TODO
def Histogram(emin, emax, nbin, A, B, C, subplts=0, lblA='', lblB='', lblC='', lbl='', Title='', totalplt=1, sm_flag=1, pltflag=1, nu=4, edge='none', edgeA='none', edgeB='none', edgeC='none', writeflag=0, flname=''):
    # find a better way to create bins, so that 0 meV is always a bin boundary
    bin_len = 3 #TODO
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
            sA = np.zeros([len(hTot[0])])
            sB = np.zeros([len(hTot[0])])
            sC = np.zeros([len(hTot[0])])
            for k in range(0,len(hTot[0])):
                sA[k] = smooth.GaussSmoothing(len(hA[0]), k, hA[0], edge=edgeA, nu=nu, X=hTot[1][:-1])
                sB[k] = smooth.GaussSmoothing(len(hB[0]), k, hB[0], edge=edgeB, nu=nu, X=hTot[1][:-1])
                sC[k] = smooth.GaussSmoothing(len(hC[0]), k, hC[0], edge=edgeC, nu=nu, X=hTot[1][:-1])
            if(pltflag == 1):
                plt.plot(hA[1][:-1], 1e1*(sA) / norm, label=lblA)
                plt.plot(hB[1][:-1], 1e1*(sB) / norm, label=lblB)
                plt.plot(hC[1][:-1], 1e1*(sC) / norm ,label=lblC)
        else:
            if(pltflag == 1):
                plt.plot(hA[1][:-1], 1e1*(hA[0] / norm), label=lblA)
                plt.plot(hB[1][:-1], 1e1*(hB[0] / norm), label=lblB)
                plt.plot(hC[1][:-1], 1e1*(hC[0] / norm), label=lblC)
    if totalplt == 1:
        if sm_flag == 1:
            sTot = np.zeros([len(hTot[0])])
            for k in range(0, len(hTot[0])):
                sTot[k] = smooth.GaussSmoothing(len(hTot[0]), k, hTot[0], edge=edge, nu=nu, X=hTot[1][:-1])
            if pltflag == 1:
                plt.plot(hTot[1][:-1], 1e1*sTot / norm, label=lbl)

        else:
            if pltflag == 1:
                plt.plot(hTot[1][:-1], 1e1*(hTot[0] / norm), label=lbl)
    if pltflag == 1:
        plt.title(Title)
        plt.legend()
        plt.grid(True)
        plt.show()
        plt.close()

    if writeflag == 1 and sm_flag == 1 and subplts == 1:
        fl = open(flname,'w')
        fl.write("#X, total, T, Q, C\n")
        for i in range(0, len(hTot[0])):
            fl.write("%e %e %e %e %e\n" %(hTot[1][i], sTot[i] / norm, sA[i] / norm, sB[i] / norm, sC[i] / norm))
        fl.close()
    elif writeflag == 1:
        print('Can only write to file with activated subplots and smoothing')


    if subplts == 1:
        if sm_flag ==1:
            return hTot[1], sTot*1e1, norm, sA*1e1, sB*1e1, sC*1e1
        else:
            return hTot[1], hTot[0]*1e1, norm, hA[0]*1e1, hB[0]*1e1, hC[0]*1e1
    else:
        if sm_flag == 1:
            return hTot[1], sTot*1e1, norm, [], [], []
        else:
            return hTot[1], hTot[0]*1e1, norm, [], [], []


def EnergyDistribution():
    return 0
