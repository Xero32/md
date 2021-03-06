####################################################################################################################
# Evaluate MD data
#
#   Track particle trajectories, calculate transition rates, compute analytical
#   solution to the Rate Equation Model
#
#   Calculate Energy Distribution, Angular Distribution
#
#
#   How To: python3 mdeval.py <Angle> <Temperature> <[--Options]>
#                       Options:    --hlrn
#                                   --nrg (float)
#                                   --start (int) [set start time for averaging of Transition Rates]
#                                   --end (int) [set end time for averaging of Transition rates]
#                                   --fp (int) [set temporal fixpoint for initial value problem in analytical solution]
#
#   Program reads files from directories    /home/<user>/lammps/111/a{angle}t{temperature}e{energy}/{job}
#   or alternatively:                       /home/<user>/lammps/111/HLRN/a{angle}t{temperature}e{energy}/{job}
#   reads '.dat' and '.lammpstrj' files, where the file names themselves are
#   the number of trajectory in the corresponding job directory, in numerical order, beginning at '1.dat'
#
#   Imports the following modules:
#       'cfg':      #NB# configuration file, where global parameters are set!
#                   e.g. save/write/plot flags for various plotting functions
#                   parameter specific variables are initialized here, but set in 'params'
#       'plot':     gives specific functions for plotting populations, transition populations, transition rates,
#                   angular distribution, energy histograms/distribution, analytical solution
#       'stats':    scaling/normalization of particle populations, various averaging/mean and standard deviation calculations,
#                   calculation of analytical solutions (eigenvalues, coefficients) and residence time, as well as
#                   initial sticking, TMAC
#       'bounce':   handles bounce events: counting and analyzing bounces, which includes evaluating particle state,
#                   transitions, transition rates, reflection angle, bounce-time resolution
#       'smooth':   includes gaussian smoothing, and other procedures (savitzky-golay filtering of trans. populations is done in 'plot')
#       'nrg':      energy distributions (normal, parallel, total, kinetic, potential, energy loss)
#       'params':   processes input arguments and parameter specific configurations, file reading functions are also found here
#       'main_mdtraj.py': main executable file, gives overview of sequenttial execution
#mdeval.py
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy, math, sys, time, argparse
import plot, bounce, nrg, smooth, stats, params
import cfg as g
import matplotlib
from pathlib import Path
#matplotlib.rcParams.update({'font.size': 16})
def printStd(a,b):
    print("%.4f +/- %.4f" %(a,b))

###### Parse Arguments
startps = 0
endps = 0
start = 0
end = 0
fp = 0
parser = argparse.ArgumentParser()
startps, endps, angle, temperature, energy, g.HLRN, fp = params.InpParams(parser)
###### Constants
kB = 1.3806503e-23
e0 = 1.60217662e-19
pi = 3.14159265359
###### Possible Input Parameters
Angle  = {'0':'0.00', '30':'0.52', '45':'0.80', '60':'1.05'}
Temperature={'80':'80', '120':'120', '160':'160', '190':'190', '220':'220', '240':'240', '270':'270', '300':'300'}
Energy = {'0':'12.95', '30':'16.027', '45':'21.71', '60':'36.55', '70':'70.0', '25.90':'25.90', '51.80':'51.80', '103.60':'103.60'}
te = 15 // 0.025 #default equilibration time
t_eq = int(fp * 0.025)+1

###### Chosen Input Parameters
params.Parameters(temperature, energy)
name = 'a'+Angle[angle]+'t'+Temperature[temperature]+'e'+Energy[energy]         # shorthand notation "a<rad>t<K>e<[float]meV>"
nameS = 'a'+angle+'t'+Temperature[temperature]+'e'+Energy[energy].split('.')[0] # shorthand "a<deg>t<K>e<[int]meV>"
nameAng = angle + u"\u00b0, " + temperature + ' K, ' + Energy[energy] + ' meV'  # notation for titles "[deg], [temp] K, [energy] meV"
''' # Ar-Pt analysis
g.jobs = (7331694, 7331705, 7331706, 7331713, 7331723, 7344079, 7344080, 7344081)
name = 'Pta'+Angle[angle]+'t'+Temperature[temperature]+'e'+Energy[energy]
'''
surface = '111'
in_folder = g.HOME + '/lammps/' + surface + '/' + name + '/'
from pathlib import Path
home = str(Path.home())
in_folder_hlrn = g.HOME + '/lammps/' + surface + '/HLRN/' + name + '/'
Temp = int(temperature)

###### Declarations
DT = 0.00025                                                                    # simulation timestep
dt = DT * 100                                                                   # data is written out only every 100 steps
maxps = -1000                                                                   # initialize the maximal simulated time in picoseconds
g.NUM_OF_JOBS = len(g.jobs)                                                     # we divide trajectory files onto multiple jobs, this is the job count
#SCALING = 1.0/(NUM_OF_JOBS*NUM_OF_TRAJ)
SCALING = 1
num = 200                                                                       # Compression Bins
num2 = 50                                                                       # 2nd Compression Bins (for Transition Rates)
wl = 0                                                                          # Compression Window Length
wl2 = 0                                                                         # 2nd Compression Window Length
Time = np.arange(0,g.TIMESTEPS_HLRN,1)                                          # time axis
State = np.zeros([3,g.TIMESTEPS_HLRN])                                          # T, Q, C states
StateComp = np.zeros([3, num])                                                  # T, Q, C states binned/compressed
Trans = np.zeros([6,g.TIMESTEPS_HLRN])                                          # QT, CT, TQ, CQ, TC, QC transition populations
TRate = np.zeros([6,num])                                                       # QT, CT, TQ, CQ, TC, QC trasnition rates
Refl = np.zeros([g.TIMESTEPS_HLRN])                                             # number of reflected particles
Bounce = np.zeros([13,80,g.NUM_OF_JOBS*g.NUM_OF_TRAJ])                          # contains information about state, time, energy of each collision ("bounce")
# nb, tm, s1, t1, E_xy.i, E_z.i, V.i, s2, t2, E_xy.f, E_z.f, V.f, refl_flag
R = np.zeros([2,2,3])                                                           # Transition Matrix
Rstd = np.zeros([2,2,3])                                                        # Standard Deviation of Transition Matrix
# energy arrays
NRG_Array = np.zeros([6, g.nbnc+1, 2, 500])                                     # contains all energy information
NormArr = np.zeros([6,g.nbnc+1])                                                # contains sclaing factor for normalizing energy distributions
delta_T = np.zeros([g.NUM_OF_JOBS*g.NUM_OF_TRAJ])                               # state dependent energy content, delta E, T state
delta_Q = np.zeros([g.NUM_OF_JOBS*g.NUM_OF_TRAJ])                               # state dependent energy content, delta E, Q state
delta_C = np.zeros([g.NUM_OF_JOBS*g.NUM_OF_TRAJ])                               # state dependent energy content, delta E, C state
E_T = np.zeros([5, g.NUM_OF_JOBS*g.NUM_OF_TRAJ])                                # state dependent energy content, total E, T state
E_Q = np.zeros([5, g.NUM_OF_JOBS*g.NUM_OF_TRAJ])                                # state dependent energy content, total E, Q state
E_C = np.zeros([5, g.NUM_OF_JOBS*g.NUM_OF_TRAJ])                                # state dependent energy content, total E, C state
paraAvg = np.zeros([2,100])                                                     # contains distributions' average value and standard deviation for parallel energy
normAvg = np.zeros([2,100])                                                     # contains distributions' average value and standard deviation for norm energy
kinAvg  = np.zeros([2,100])                                                     # contains distributions' average value and standard deviation for kinetic energy
deltaAvg  = np.zeros([2,100])                                                   # contains distributions' average value and standard deviation for exchange energy
binsArr = [141,80,100,120,145,120]  # deltaE, paraE, normE, potE, totE, kinE    # how many bins for which histogram/distribution
# misc
hlrnflag = 0
globalcounter = 0
falsecounter = 0
TOTAL = g.NUM_OF_JOBS*g.NUM_OF_TRAJ
bounceavg = np.zeros([TOTAL])
horizontal = np.zeros([4,len(Time)])


def ReadfileFn():
    global name
    global globalcounter
    global falsecounter
    global TOTAL
    global maxps
    global dt

    print('Evaluating Data for', name)
    STARTTIME = time.time()
    ###### Read Trajectory Files
    for ctr, jb in enumerate(g.jobs):
        if jb == 0: continue
        for d in range(1, g.NUM_OF_TRAJ+1):
            # read files
            Traj, g.TIMESTEPS, hlrnflag = params.Readfiles(ctr, d, jb, name, in_folder, in_folder_hlrn) #, TIMESTEPS_HLRN, TIMESTEPS_RZ, NUM_OF_TRAJ_RZ)
            maxps = g.TIMESTEPS * dt

            g.TRAJ_COUNTER += 1
            if hlrnflag == 1:
                g.HLRN_COUNTER += 1
            jj = globalcounter
            globalcounter += 1
            # count and analyze bounces
            nb = bounce.Countbounces(Traj[3,:], Traj[1,:], Traj[0,:], Traj[2,:], Bounce[:,:,jj], jj, g.TIMESTEPS)

            # sort out erratic trajectories
            if nb == -9999:
                falsecounter += 1
                TOTAL -= 1
                g.TRAJ_COUNTER -= 1
                if hlrnflag == 1:
                    g.HLRN_COUNTER -= 1
                continue

            # create state populations and transition populations
            for i in range(0,nb):
                try:
                    flag = bounce.Population(i, Bounce[:,i-1,jj], Bounce[:,i,jj], Bounce[:,i+1,jj], State, Trans, Refl, 1., jj, g.TIMESTEPS)
                except:
                    flag = bounce.Population(i, Bounce[:,i,jj], Bounce[:,i,jj], Bounce[:,i+1,jj], State, Trans, Refl, 1., jj, g.TIMESTEPS)
                # look for errors
                if flag == -6000 or flag == -5000:
                    TOTAL -= 1
                    g.TRAJ_COUNTER -= 1
                    falsecounter += 1
                    if hlrnflag == 1:
                        g.HLRN_COUNTER -= 1
                        falsecounter += 1
                    bounceavg[jj] = 10000
                    print(str(flag))
                    break
                elif flag == -3000:
                    bounceavg[jj] = nb
                    break
                if flag == 1:
                    Bounce[12,i,jj] = flag
                    break
            bounceavg[jj] = nb
            ti = int(Bounce[3,0,jj])
            tf = int(Bounce[8,nb-1,jj])
            if tf >= g.TIMESTEPS:
                tf -= 1
            if Bounce[7,nb-1,jj] == 1:
                pass
                # TODO calc TMAC or angular distributions


        #end for (trajectory)
    #end for (job list)
    print("Found %d miscounted trajectories in %d total trajectories and deleted them." %(falsecounter,g.TRAJ_COUNTER+falsecounter) )
    print("TIMESTEPS: ", g.TIMESTEPS)
    print("HLRN trajectories: ", g.HLRN_COUNTER)

    # rescale populations so that we are inside [0.0, 1.0]
    stats.Scaling(State, Trans, Refl)
    bncCumu = bounce.bncHist(bounceavg)

    ENDTTIME = time.time()
    #print("The loops took %f seconds" %(ENDTTIME-STARTTIME))



###### TMAC
# avgtmac = stats.median(Tmac, fltr=1, fltrval=0)
# stdtmac = stats.TMACstd(Tmac, avgtmac)
#printStd(avgtmac, stdtmac)
'''
fn = open("TMAC.dat", "a")
fn.write("%d %d %f %f %f\n" %(int(angle), int(temperature), float(Energy[energy]), avgtmac, stdtmac))
fn.close()
'''

"""
###### Angular Distribution
NB = [1,2,4,6,8,10,12,14]
assert(len(NB) == 8)
lblNB = ['0','2', '4', '6', '8', '10', '12', '14']
nbin = 36
theta0 /= thetactr
# ThetaEQ = bounce.thetaEQ(Theta, NB, g.TRAJ_COUNTER+falsecounter)

locs = [0, pi/12, pi/6, pi/4, pi/3, 5*pi/12, pi/2]
labels = ['0', r'$\frac{\pi}{12}$', r'$\frac{\pi}{6}$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{3}$', r'$\frac{5\pi}{12}$', r'$\frac{\pi}{2}$']
#labels = ['0°', '15°', '30°', '45°', '60°', '75°', '90']
if g.P_ANG != 0:
    plt.xticks(locs, labels)
if g.HLRN == 1:
    svnameAng = g.HOME+'/lammps/newplot/ang/'+nameS+'HLRNAngDistr.pdf'
else:
    svnameAng = g.HOME+'/lammps/newplot/ang/'+nameS+'AngDistr.pdf'
"""
'''
plot.Distr(0.0, np.pi/2 + np.pi/(2*nbin), nbin, ThetaEQ, subplts=0, lbl=['1','2', '4', '6', '8', '10', '12', '14'], numTraj=g.TRAJ_COUNTER+falsecounter,  NB=NB,
 Title='Angular Distribution, ' + nameAng, binarr=[], avg=0, std=0, xlbl='$\Theta$ / $rad$', ylbl='$p(\Theta)$ / rad$^{-1}$', a=angle, t=temperature, e=Energy[energy],
 saveflag=g.S_ANG, savename=svnameAng, subscale=1, totalplt=1, sm_flag=g.G_ANG, pltflag=g.P_ANG, nu=2, edge='none', writeflag=g.W_ANG, flname=g.HOME+'/lammps/AngDistr.dat', nb=0, action=-1, scl=1, vertval=[], vertls=[])
'''
###### Calculate Sticking Probability
def StickingFn():
    tD = 10 // 0.025
    StickCoeff = stats.InitSticking(State[0], State[1], tD, int(angle),int(temperature), float(Energy[energy]), g.TRAJ_COUNTER, writeflag=g.W_STICK, filename='DelayStickingProb.dat', printflag=1)


###### Particle Populations
def PlotPopulations():
    plt.clf()
    plt.cla()

    plot.Populations(angle, temperature, Energy[energy], Time, State[0], State[1], State[2], smflag=g.G_POP, pltflag=g.P_POP, hlrn=g.HLRN, Title=nameAng)
###### Compress Population for better statistics and easier handling
def Binning():
    global num
    global wl

    for i in range(0,3):
        StateComp[i], wl, num = smooth.Compress(State[i], num=num)

    for i in range(0,6):
        Trans[i, :num], wl, num = smooth.Compress(Trans[i], num=num)

    TimePrime = np.arange(0,g.TIMESTEPS_HLRN, wl) # Compressed Timescale
    return TimePrime

def PlotTransitionPopulations(TimePrime):
    if int(temperature) < 160:
        nu = 25 # Savitzky-Golay Smoothing Parameter
    else:
        nu = 45
    plot.TransitionPopulations(angle, temperature, Energy[energy],
    TimePrime, Trans[0, :num], Trans[2, :num], Trans[1, :num], Trans[3, :num], smflag=g.G_TRANS, pltflag=g.P_TRANS, hlrn=g.HLRN, Title=nameAng, nu=nu)

def Calc_BinTransitionRates():
    global num2
    global wl
    global wl2
    global maxps

    print("Calculate Transition Rates")
    d = 0
    startdiff = 5
    for i in range(0,6):
        TRate[i] = bounce.DifferentiateT(1, StateComp[d], Trans[i, :num], TRate[i], wl, startdiff, maxps)
        d += i % 2

    ###### Further Compress Transition Rates
    for i in range(0,6):
        TRate[i, :num2], wl2, num2 = smooth.Compress(TRate[i], num=num2)

    TimePrime2 = np.arange(0, g.TIMESTEPS_HLRN, wl*wl2) # 2nd Compressed Timescale
    return TimePrime2

###### Plot Transition Rates, calculate avg value
def PlotTransitionRates(TimePrime2):
    global maxps, startps, endps
    global num2
    global start, end, fp
    start = int(startps * num2 / maxps)
    end = int(endps * num2 / maxps)
    avgflag = g.P_T_AVG
    if int(temperature) < 160:
        # avgflag = 0
        nu = 4
    else:
        # avgflag = 1
        nu = 8
    plot.TransitionRate(angle, temperature, Energy[energy], TimePrime2, TRate[0, :num2], TRate[2, :num2], TRate[1, :num2], TRate[3, :num2], lblA=r'$T_{QT}$', lblB=r'$T_{TQ}$', lblC=r'$T_{CT}$', lblD=r'$T_{CQ}$',
    smflag=1, pltflag=g.P_T, ylbl='Transition Rate / t\u2080\u207B\u00B9', avgflag=avgflag, start=start, end=end, hlrn=g.HLRN, Title=nameAng)


###### Setup Analytical Solution
def AnalyticalSolution(writeflag=g.W_PARAM):
    global maxps, num2
    global start, end, fp
    global Angle, Temp, Energy, angle, energy, surface
    print("maxps: ", maxps)
    print("start_index: ", start, " end_index: ", end)
    multiplier=3 # Sets how long simulation time is to be exceeded by analytical solution
    N, l1, l2, delL, c1, dc1, c2, dc2, Time2 = stats.CalcSolution(fp, State, start, end, R, Rstd, TRate, num2, multiplier=multiplier)
    ###### Calc Residence Time
    tR = stats.ResTime(l1, N[0][fp], N[1][fp], R[:,:,0])
    tRstd = stats.StdResTime(l1, delL, N[0][fp], N[1][fp], R[:,:,0])
    print('\ntR:')
    printStd(tR, tRstd)
    print('\n')
    avgt0 = stats.AvgRate(start, end, TRate[0])
    stdt0 = stats.StdRate(start, end, TRate[0], avgt0)
    avgt1 = stats.AvgRate(start, end, TRate[1])
    stdt1 = stats.StdRate(start, end, TRate[1], avgt1)
    avgt2 = stats.AvgRate(start, end, TRate[2])
    stdt2 = stats.StdRate(start, end, TRate[2], avgt2)
    avgt3 = stats.AvgRate(start, end, TRate[3])
    stdt3 = stats.StdRate(start, end, TRate[3], avgt3)

    if (writeflag != 0):
        ## write eigenvalues, transition matrix and coefficients of analytical solution to extra file
        ## for further analysis for particle fluxes
        home = str(Path.home())
        paramname = str("Single_a%st%de%s.dat" % (Angle[angle], Temp, Energy[energy]))
        fname = home + "/lammps/" + surface + "/" + paramname
        f = open(fname, 'w')
        f.write("# Eigenvalues and Coefficients for Single Particle Solution\n")
        f.write("lambda1 = %.16f\n" % l1)
        f.write("lambda2 = %.16f\n" % l2)
        f.write("c1 = %.16f\n" % c1)
        f.write("c2 = %.16f\n" % c2)
        f.write("# Elements of Transition Matrix\n")
        f.write("R11 = %.16f\n" % R[0,0,0])
        f.write("R12 = %.16f\n" % R[0,1,0])
        f.write("R21 = %.16f\n" % R[1,0,0])
        f.write("R22 = %.16f\n" % R[1,1,0])
        f.write("N1 = %.16f\n" % N[0][fp])
        f.write("N2 = %.16f\n" % N[1][fp])
        f.write("# Equilibration time in ps\n")
        f.write("te = %d\n" % t_eq)
        f.write('# Transition rates [t_0^{-1}]\n')
        f.write('T_QT = %.16f\n' % avgt0)
        f.write('T_CT = %.16f\n' % avgt1)
        f.write('T_TQ = %.16f\n' % avgt2)
        f.write('T_CQ = %.16f\n' % avgt3)
        f.close()


    ###### Give out Average Transition Rates for further Analysis
    horizontal[0] = np.array([R[1,0,0] for i in range(len(Time))]) #T_QT
    horizontal[1] = np.array([R[0,1,0] for i in range(len(Time))]) #T_TQ
    horizontal[2] = np.array([-R[0,0,0]-R[1,0,0] for i in range(len(Time))]) #T_CT
    horizontal[3] = np.array([-R[1,1,0]-R[0,1,0] for i in range(len(Time))]) #T_CQ


    return N, Time2, multiplier

def Write_Const_TransitionRates():
    if g.W_T != 0:
        f = open('TransitionRate.txt','a')
        f.write("%d %d %f, %f %f %f %f %f %f %f %f\n" % (int(angle), int(Temp), float(Energy[energy]), avgt0, stdt0, avgt1, stdt1, avgt2, stdt2, avgt3, stdt3))
        f.close()

###### Equilibration Condition
def Check_EquilibrationCondition(flag):
    if flag == True:
        left = StateComp[0,-1] * horizontal[0,-1]
        right = StateComp[1,-1] * horizontal[2,-1] + StateComp[1,-1] * horizontal[3,-1]
        print('N_T * T_QT \t\tvs.\t N_Q * (T_TQ + T_CQ)')
        print(left, '\tvs.\t', right)

'''
fl = open("ResTime.dat", 'a')
fl.write("%2d %3d %.3f %.4f %.4f %.4f %.4f\n" %(int(angle), int(temperature), float(Energy[energy]), tR, tRstd, tRsimple/tR, tP))
fl.close()
'''
def PlotSolution(N, Time2, TimePrime, multiplier, te=0):
    global dt
    global t_eq, num, maxps
    te_index = int((t_eq / maxps) * num)
    if g.HLRN == 1:
        nameA = 'a'+angle+'t'+Temperature[temperature]+'e'+Energy[energy]+'HLRN'
    else:
        nameA = 'a'+angle+'t'+Temperature[temperature]+'e'+Energy[energy]
    home = str(Path.home())
    path = home + "/lammps/111/"

    plot.WritePlot(X=TimePrime[:te_index]*0.025, Y=StateComp[0,:te_index], name=path+nameA+"InitTrapped", xlabel='# t / ps', ylabel='Population fraction', saveflag=True, header=True, action='w')
    plot.WritePlot(X=TimePrime[:te_index]*0.025, Y=StateComp[1,:te_index], name=path+nameA+"InitQuasiTr", xlabel='# t / ps', ylabel='Population fraction', saveflag=True, header=True, action='w')
    plot.WritePlot(X=TimePrime[:te_index]*0.025, Y=StateComp[2,:te_index], name=path+nameA+"InitCont", xlabel='# t / ps', ylabel='Population fraction', saveflag=True, header=True, action='w')
    plot.Solution(angle, temperature, Energy[energy], N, Time2, TimePrime, StateComp, Title=nameAng, pltflag=g.P_SOL, maxtime=int(g.TIMESTEPS*multiplier)*dt)

def EnergyHandling():
    # arrays for mean energy value, needed for equilibration analysis
    wf = g.W_NRG #writeflag
    pf = 0

    if g.HLRN == 1:
        nameS1 = nameS+'HLRN'
    else:
        nameS1 = nameS

    filename = g.HOME+"/lammps/111/nrg/" + nameS1 + 'E.dat'

    binsArrL =[]
    for i in range(6):
        binsArrL.append(binsArr[i]-1)

    print("Evaluate Energy Distributions")
    for n in range(0,g.nbnc+1,1):
        #sprint(str(n))
        #TODO
        #filename = in_folder + 'deltaE' + str(n) + '.dat'
        act = 0
        t,q,c = nrg.EnergyLoss(n, delta_T, delta_Q, delta_C, g.NUM_OF_JOBS*g.NUM_OF_TRAJ, Bounce[4:7,:,:], Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
        delE, distrDelE, NormArr[act,n], a1,a2,a3, average, stdev = plot.Histogram(-130,150,binsArr[act], delta_T[:t], delta_Q[:q], delta_C[:c], subplts=1, lblA='delta T', lblB='delta Q', lblC='delta C', lbl='Total delta E',
        pltflag=pf, sm_flag=1, nu=2, writeflag=wf, flname=filename, nb=n, action=act, scl=1e1, avg=1, std=1)
        #nrg.PlotHist(-120,140,80, delta_T[:t], delta_Q[:q], delta_C[:c], subplts=1, lblA='delta T', lblB='delta Q', lblC='delta C', lbl='Total delta E', pltflag=1, sm_flag=1)
        deltaAvg[0,n] = average
        deltaAvg[1,n] = stdev

        NRG_Array[act, n, 0, :binsArr[act]] = delE              # x-axis
        NRG_Array[act, n, 1, :binsArr[act]-1] = distrDelE       # y-axis

        #filename = in_folder + 'paraE' + str(n) + '.dat'
        act += 1
        ## normal bin setting: 0,80,80: change last number (bin number) in binsArr
        t,q,c = nrg.paraEnergyDistr(n, E_T[0,:], E_Q[0,:], E_C[0,:], g.NUM_OF_JOBS*g.NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
        paraE, distrParaE, NormArr[act,n], a1, a2, a3, average, stdev = plot.Histogram(0,80,binsArr[act], E_T[0][:t], E_Q[0][:q], E_C[0][:c], lblA='T', lblB='Q', lblC='C', lbl='Parallel Energy Distr',
        edge='pos', edgeA='pos', edgeB='pos', edgeC='pos',
        subplts=1, sm_flag=1, pltflag=pf, nu=1+math.sin(n/30*math.pi/2), writeflag=wf, flname=filename, avg=1, std=1, nb=n, action=act, scl=1e1)
        # with the sine function we increase the smoothing with growing n, as then the number of existing trajectories decreases
        paraAvg[0,n] = average

        NRG_Array[act,n,0,:binsArr[act]] = paraE
        NRG_Array[act,n,1,:binsArr[act]-1] = distrParaE

        #filename = in_folder + 'normE' + str(n) + '.dat'
        act += 1
        ## normal bin setting: 0,150,100 maybe
        t,q,c = nrg.normEnergyDistr(n, E_T[1,:], E_Q[1,:], E_C[1,:], g.NUM_OF_JOBS*g.NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
        normE, distrNormE, NormArr[act,n], a1, a2, a3, average, stdev = plot.Histogram(0,150,binsArr[act], E_T[1][:t], E_Q[1][:q], E_C[1][:c], subplts=1, lblA='T', lblB='Q', lblC='C', lbl='Normal Energy Distr',
        edge='pos', edgeA='pos', edgeB='pos', edgeC='pos', avg=1, std=1,
        sm_flag=1, pltflag=pf, nu=2, writeflag=wf, flname=filename, nb=n, action=act, scl=1e1)
        normAvg[0,n] = average

        NRG_Array[act,n,0,:binsArr[act]] = normE
        NRG_Array[act,n,1,:binsArr[act]-1] = distrNormE

        #filename = in_folder + 'potE' + str(n) + '.dat'
        act += 1
        t,q,c = nrg.potEnergyDistr(n, E_T[2,:], E_Q[2,:], E_C[2,:], g.NUM_OF_JOBS*g.NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
        potE, distrPotE, NormArr[act,n], *waste = plot.Histogram(-140,-60,binsArr[act], E_T[2][:t], E_Q[2][:q], E_C[2][:c], subplts=1, lblA='T', lblB='Q', lblC='C', lbl='Potential Energy Distr', sm_flag=1, pltflag=pf, nu=2, writeflag=wf, flname=filename,
        nb=n, action=act, scl=1e1)
        NRG_Array[act,n,0,:binsArr[act]] = potE
        NRG_Array[act,n,1,:binsArr[act]-1] = distrPotE

        #filename = in_folder + 'totE' + str(n) + '.dat'
        act += 1
        t,q,c = nrg.totEnergyDistr(n, E_T[3,:], E_Q[3,:], E_C[3,:], g.NUM_OF_JOBS*g.NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
        totE, distrTotE, NormArr[act,n], *waste = plot.Histogram(-180,180,binsArr[act], E_T[3][:t], E_Q[3][:q], E_C[3][:c], subplts=1, lblA='T', lblB='Q', lblC='C', lbl='Total Energy Distr', edge='none', edgeA='neg', edgeB='pos', edgeC='pos',
        sm_flag=1, pltflag=pf, nu=2, writeflag=wf, flname=filename, nb=n, action=act, scl=1e1)

        NRG_Array[act,n,0,:binsArr[act]] = totE
        NRG_Array[act,n,1,:binsArr[act]-1] = distrTotE

        act += 1
        t,q,c = nrg.kinEnergyDistr(n, E_T[4,:], E_Q[4,:], E_C[4,:], g.NUM_OF_JOBS*g.NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
        kinE, distrKinE, NormArr[act,n], a1,a2,a3, average, stdev = plot.Histogram(0,200,binsArr[act], E_T[4][:t], E_Q[4][:q], E_C[4][:c], lblA='T', lblB='Q', lblC='C', lbl='Kinetic Energy Distr',
        edge='pos', edgeA='pos', edgeB='pos', edgeC='pos', Title='KE, nb= '+str(n),
        subplts=1, sm_flag=1, pltflag=pf, nu=2, writeflag=wf, flname=filename, avg=1, std=1, nb=n, action=act, scl=1e1)
        kinAvg[0,n] = average
        NRG_Array[act,n,0,:binsArr[act]] = kinE
        NRG_Array[act,n,1,:binsArr[act]-1] = distrKinE

    return nameS1


###### plot energy distributions
def PlotEnergyDistr(nameS1):
    act = 0
    Bounces, axes = plot.ChooseParams(Temp, int(angle))
    Action = ['deltaE', 'paraE', 'normE', 'potE', 'totE', 'kinE']
    Title = ['Energy Loss Distribution', 'Parallel Energy Distribution', 'Normal Energy Distribution', 'Potential Energy Distribution', 'Total Energy Distribution', 'Kinetic Energy Distribution']
    Xlbl = [r'$\Delta E$ / meV', r'$E_{||}$ / meV', r'$E_{\bot}$ / meV', r'$V$ / meV', r'$E$ / meV', r'$E_{kin}$ / meV']
    Ylbl = [r'P($\Delta E ,n_b$) / $10^{-1}$ meV$^{-1}$', r'P($E_{||} ,n_b$) / $10^{-1}$ meV$^{-1}$', r'P($E_{\bot} ,n_b$) / $10^{-1}$ meV$^{-1}$', r'P($V ,n_b$) / $10^{-1}$ meV$^{-1}$',
     r'P($E ,n_b$) / $10^{-1}$ meV$^{-1}$', r'P($E_{kin} ,n_b$) / $10^{-1}$ meV$^{-1}$']

    print(nameS1 + ': ' + str(deltaAvg[0,Bounces[act][-1]]) + ' +/- ' + str(deltaAvg[1,Bounces[act][-1]]))
    pflag = g.P_NRG
    for act in range(6):
        plot.distrEnergy(Temp, act, Bounces[act], NRG_Array, NormArr[act], binsArr[act]-1, xlbl=Xlbl[act], ylbl=Ylbl[act],
          ax=axes[act], saveflag=g.S_NRG, svname="/home/becker/lammps/newplot/nrg/"+Action[act]+'/'+nameS1+Action[act]+".pdf", pltflag=pflag, Title=Title[act], figname=nameS1)


def CalcTrajectoryDecay():
    #evaluate how many bounces typically appear in a trajectory with fixed timesteps at Nt = 120000
    pltmax = g.nbnc+10
    sequence = np.arange(0,pltmax+1)
    plt.show(block=False)
    bouncehist = np.histogram(bounceavg, density=False, bins=sequence)

    #look at how the number of trajectories diminishes as nb grows
    bouncecumulative = np.cumsum(bouncehist[0])
    return bouncecumulative

def Std_EnergyOverTime(bouncecumulative):
    #calc StdDev for energy development over time
    for n in range(0,g.nbnc+1):
        paraAvg[1,n] = stats.stdDev(Bounce[9,n,:], paraAvg[0,n], (g.NUM_OF_JOBS*g.NUM_OF_TRAJ)-bouncecumulative[n])
        normAvg[1,n] = stats.stdDev(Bounce[10,n,:], normAvg[0,n], (g.NUM_OF_JOBS*g.NUM_OF_TRAJ)-bouncecumulative[n])
        kinAvg[1,n] = stats.stdDev(Bounce[10,n,:]+Bounce[9,n,:], kinAvg[0,n], (g.NUM_OF_JOBS*g.NUM_OF_TRAJ)-bouncecumulative[n])

######## plot mean Energy over time here
def Plot_MeanEnergyOverTime():
    pltmax = g.nbnc+10
    sequence = np.arange(0,pltmax+1)
    tempArr = [int(temperature) * kB / e0 *1e3 for i in range(len(sequence))]


    plt.xlabel("n")
    plt.ylabel("E / meV")
    plt.errorbar(sequence, paraAvg[0,:pltmax+1], paraAvg[1,:pltmax+1], color='red', label='parallel kin. Energy / meV')
    plt.errorbar(sequence, normAvg[0,:pltmax+1], normAvg[1,:pltmax+1], color='black', label = 'normal kin. Energy / meV')
    plt.errorbar(sequence, kinAvg[0,:pltmax+1], kinAvg[1,:pltmax+1], color='green', label = 'total kin. Energy / meV')
    plt.plot(sequence, tempArr, color='blue', ls='--')
    plt.legend()
    plt.show(block=False)

def Find_EquilibrationTime():
    pltmax = g.nbnc+10
    sequence = np.arange(0,pltmax+1)
    ## fitness function, which takes into account deviations from neighboing points and surface temp, and takes
    ## into account number of bounces and errorbars
    ##
    ## aim is to find the time of maximal thermalization of the adatom
    def fct(inArr, errArr, fixbounce, temp):
        lossArr = np.full((len(inArr)), 9999.)
        minimum = 9999.
        minimumplace = -10
        for i in range(0,len(inArr)-2):
            t1 = math.fabs(inArr[i] - inArr[i+1])
            t2 = math.fabs(inArr[i] - inArr[i-1])
            t3 = errArr[i]
            t4 = math.fabs(float(fixbounce - i))
            t5 = math.fabs(float(temp) * kB / e0 * 1000. - inArr[i])
            # term6 = tRate[]
            lossArr[i] = (t1 * 1.2 + t2) * 5. + t3 * 3. + t4 * 0.1 + t5 * 0.04
        for i in range(5, g.nbnc):
            if lossArr[i+1] - lossArr[i] > 50:
                break
            if lossArr[i] < minimum:
                minimum = lossArr[i]
                minimumplace = i
        return minimumplace, lossArr

    minplace1, lossArr = fct(paraAvg[0], paraAvg[1], g.nbnc - 5, int(temperature))
    print(minplace1)
    plt.plot(sequence, lossArr[:pltmax+1], label='para energy')

    minplace2, lossArr = fct(normAvg[0], normAvg[1], g.nbnc - 5, int(temperature))
    print(minplace2)
    plt.plot(sequence, lossArr[:pltmax+1], label='norm energy')

    minplace3, lossArr = fct(kinAvg[0], kinAvg[1], g.nbnc - 5, int(temperature))
    print(minplace3)
    plt.plot(sequence, lossArr[:pltmax+1], label='total energy')
    plt.legend()
    plt.show(block=False)
    plt.clf()
    plt.cla()
    plt.close()

    effPara = paraAvg[0,minplace1] / kB * e0 /1e3
    effNorm = normAvg[0,minplace2] / kB * e0 /1e3
    effTemp = effPara + effNorm
    effNRG = paraAvg[0,minplace1] + normAvg[0,minplace2]
    print(effTemp)
    print(effNRG)

def Write_ResidenceTime():
    flagResTime = 0
    if flagResTime != 0:
        fname = "ResidenceTimeNorm.dat"
        f = open(fname, 'a')
        f.write("%d %d %f %f %f %f\n" %(int(angle), int(temperature), effNorm, float(Energy[energy]), tR, tRstd))
        f.close()


def PlotTrajectoryDecay(bouncecumulative):
    global falsecounter
    pltmax = g.nbnc+10
    sequence = np.arange(0,pltmax+1)
    mean = ((g.NUM_OF_JOBS*g.NUM_OF_TRAJ + falsecounter)-bouncecumulative[0]) / 2
    if g.HLRN == 1:
        name = 'a'+angle+'t'+temperature+'e'+Energy[energy].split('.')[0]+'HLRN'
    else:
        name = 'a'+angle+'t'+temperature+'e'+Energy[energy].split('.')[0]
    svname = g.HOME+'/lammps/newplot/TrajDecr/' + name + 'TrajDecr.pdf'
    if g.P_DECR != 0:
        plt.plot(np.arange(0,pltmax), ((g.NUM_OF_JOBS*g.NUM_OF_TRAJ)-bouncecumulative), color='k')
        # plt.plot(np.arange(0,pltmax), [mean for i in range(pltmax)], 'k:')
        plot.SetupPlot('Number of Particles, ' + nameAng,  r'$n_b$', 'Number of Particles', grid=True, Block=1-g.S_DECR, saveflag=g.S_DECR, savename=svname, sz=12, legend=0)

def CalcAvgBounceTime():
    global TOTAL
    #calc average bounce time + stddev
    nbnc = g.nbnc
    TimeArr = np.zeros([3,80])
    TimeArr = bounce.BncTime(Bounce, TOTAL, nbnc, g.TIMESTEPS, TimeArr)
    timebetween = bounce.TimeBetwBounces(TimeArr[:,:g.nbnc]) * 2.5e-2
    print("Time between bounces:",timebetween,'ps')
    if g.HLRN == 1:
        name = 'a'+angle+'t'+temperature+'e'+Energy[energy].split('.')[0]+'HLRN'
    else:
        name = 'a'+angle+'t'+temperature+'e'+Energy[energy].split('.')[0]
    svname = g.HOME+'/lammps/newplot/TimeBetw/' + name + 'TimeBetw.pdf'
    if g.P_TBETW != 0:
        plt.cla()
        plt.clf()
        plt.errorbar(TimeArr[0][:g.nbnc]*2.5e-2, np.arange(0,g.nbnc,1), xerr=TimeArr[1][:g.nbnc]*2.5e-2, color='k')
        plot.SetupPlot('Time of Bounces, ' + nameAng, 't / ps', r'$n_b$', grid=True, Block=1-g.S_TBETW, saveflag=g.S_TBETW, savename=svname, sz=12, legend=0)
