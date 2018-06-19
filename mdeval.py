#mdeval.py
import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import sys
import plot
import bounce
import nrg
import smooth
import stats
import argparse
import params
startps = 0
endps = 0
#TIMESTEPS, TIMESTEPS_RZ, TIMESTEPS_HLRN, TRAJ_COUNTER, HLRN_COUNTER, HLRN, NUM_OF_JOBS, NUM_OF_TRAJ, NUM_OF_TRAJ_RZ
G_LIST = np.zeros(20)
parser = argparse.ArgumentParser()
startps, endps, angle, temperature, energy, HLRN = params.InpParams(parser)

kB = 1.3806503e-23
e0 = 1.60217662e-19

Angle  = {'0':'0.00', '30':'0.52', '45':'0.80', '60':'1.05'}
Temperature={'80':'80', '120':'120', '160':'160', '190':'190', '220':'220', '240':'240', '270':'270', '300':'300'}
Energy = {'0':'12.95', '30':'16.027', '45':'21.71', '60':'36.55', '70':'70.0'}
te = 15 // 0.025 #default
TIMESTEPS_RZ = 120000 // 100
TIMESTEPS_HLRN = 120000 // 100
NUM_OF_TRAJ_RZ = 200
NUM_OF_TRAJ = 200

jobs = (7227490, 7227494, 7227495, 7227496, 7227497)
if energy != '70':
    TIMESTEPS_RZ, TIMESTEPS_HLRN, NUM_OF_TRAJ_RZ, NUM_OF_TRAJ, nbnc, jobs = params.Parameters(temperature, HLRN)

name = 'a'+Angle[angle]+'t'+Temperature[temperature]+'e'+Energy[energy]
''' Ar-Pt analysis
jobs = (7331694, 7331705, 7331706, 7331713, 7331723, 7344079, 7344080, 7344081)
name = 'Pta'+Angle[angle]+'t'+Temperature[temperature]+'e'+Energy[energy]
'''
in_folder = '/home/becker/lammps/111/' + name + '/'
#in_folder = '/home/becker/lammps/111/AuLong/' + name + '/'
Temp = int(temperature)


DT = 0.00025
dt = DT * 100
NUM_OF_JOBS = len(jobs)
#SCALING = 1.0/(NUM_OF_JOBS*NUM_OF_TRAJ)
SCALING = 1
dataflag = 0
#TOTAL = np.zeros([TIMESTEPS_HLRN])
Time = np.arange(0,TIMESTEPS_HLRN,1)
State = np.zeros([3,TIMESTEPS_HLRN]) # T, Q, C
Trans = np.zeros([6,TIMESTEPS_HLRN]) # QT, CT, TQ, CQ, TC, QC
TRate = np.zeros([6,TIMESTEPS_HLRN]) # QT, CT, TQ, CQ, TC, QC
Refl = np.zeros([TIMESTEPS_HLRN])
Bounce = np.zeros([13,80,NUM_OF_JOBS*NUM_OF_TRAJ])
# nb, tm, s1, t1, E_xy.i, E_z.i, V.i, s2, t2, E_xy.f, E_z.f, V.f, refl_flag
TRAJ_COUNTER = 0
HLRN_COUNTER = 0
hlrnflag = 0
globalcounter = 0
falsecounter = 0
TOTAL = NUM_OF_JOBS*NUM_OF_TRAJ
bounceavg = np.zeros([TOTAL])

print('Evaluating Data for', name)
for ctr, jb in enumerate(jobs):
    for d in range(1, NUM_OF_TRAJ+1):
        # read files
        Traj, TIMESTEPS, hlrnflag = params.Readfiles(ctr, d, jb, name, in_folder, TIMESTEPS_HLRN, TIMESTEPS_RZ, NUM_OF_TRAJ_RZ)

        TRAJ_COUNTER += 1
        if hlrnflag == 1:
            HLRN_COUNTER += 1
        jj = globalcounter
        globalcounter += 1

        nb = bounce.Countbounces(Traj[3,:], Traj[1,:], Traj[0,:], Traj[2,:], Bounce[:,:,jj], jj, TIMESTEPS)

        if nb == -9999:
            falsecounter += 1
            TOTAL -= 1
            TRAJ_COUNTER -= 1
            if hlrnflag == 1:
                HLRN_COUNTER -= 1
            continue

        for i in range(0,nb):
            try:
                flag = bounce.Population(i, Bounce[:,i-1,jj], Bounce[:,i,jj], Bounce[:,i+1,jj], State, Trans, Refl, 1., jj, TIMESTEPS)
            except:
                flag = bounce.Population(i, Bounce[:,i,jj], Bounce[:,i,jj], Bounce[:,i+1,jj], State, Trans, Refl, 1., jj, TIMESTEPS)
            if flag == -6000 or flag == -5000:
                TOTAL -= 1
                TRAJ_COUNTER -= 1
                falsecounter += 1
                if hlrnflag == 1:
                    HLRN_COUNTER -= 1
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

    #end for (trajectory)
#end for (job list)

print("Found %d miscounted trajectories in %d total trajectories and deleted them." %(falsecounter,TRAJ_COUNTER+falsecounter) )
print("TIMESTEPS: ", TIMESTEPS)
print("HLRN trajectories: ", HLRN_COUNTER)

stats.Scaling(State, Trans, Refl, TRAJ_COUNTER, HLRN_COUNTER, TIMESTEPS_RZ, TIMESTEPS_HLRN)

StateNew = np.zeros([3, TIMESTEPS_HLRN])
StateSum = np.zeros([3, TIMESTEPS_HLRN])
TransNew = np.zeros([6, TIMESTEPS_HLRN])
TRateE   = np.zeros([6, TIMESTEPS_HLRN])
TRateA   = np.zeros([6, TIMESTEPS_HLRN])
                #4       2              0           1
StateSum[0] = Trans[4] + Trans[2] - Trans[0] - Trans[1]
                #5            0           2           3
StateSum[1] = Trans[5] + Trans[0] - Trans[2] - Trans[3]
#C              CT          CQ          TC      QC
StateSum[2] = 1 + Trans[1] + Trans[3] - Trans[4] - Trans[5]

StickCoeff = stats.InitSticking(Refl, float(Angle[angle]),float(Temperature[temperature]), float(Energy[energy]), TRAJ_COUNTER)

plot.Populations(angle, temperature, Energy[energy], Time, State[0], State[1], State[2], smflag=1, pltflag=1, hlrn=HLRN)
num = 200
StateComp = np.zeros([3, num])
for i in range(0,3):
    StateComp[i], wl, num = smooth.Compress(State[i], num=num)

for i in range(0,6):
    Trans[i, :num], wl, num = smooth.Compress(Trans[i], num=num)

TimePrime = np.arange(0,TIMESTEPS_HLRN, wl)

plot.TransitionPopulations(angle, temperature, Energy[energy],
TimePrime, Trans[0, :num], Trans[2, :num], Trans[1, :num], Trans[3, :num], smflag=1, pltflag=1, hlrn=HLRN)


'''State = np.zeros([3,TIMESTEPS]) # T, Q, C
Trans = np.zeros([4,TIMESTEPS]) # QT, CT, TQ, CQ
TRate = np.zeros([4, TIMESTEPS]) # QT, CT, TQ, CQ, TC, QC'''
TRate = np.zeros([6, num])
print("Calculate Transition Rates")
d = 0
maxps = TIMESTEPS * dt
startdiff = 10
for i in range(0,6):
    TRate[i] = bounce.DifferentiateT(1, StateComp[d], Trans[i, :num], TRate[i], wl, startdiff, maxps)
    d += i % 2

num2 = 50
for i in range(0,6):
    TRate[i, :num2], wl2, num2 = smooth.Compress(TRate[i], num=num2)

TimePrime2 = np.arange(0, TIMESTEPS_HLRN, wl*wl2)

plot.TransitionRate(angle, temperature, Energy[energy], TimePrime2, TRate[0, :num2], TRate[2, :num2], TRate[1, :num2], TRate[3, :num2],
lblA='T_QT', lblB='T_TQ', lblC='T_CT', lblD='T_CQ', smflag=1, pltflag=0, ylbl='T / t\u2080\u207B\u00B9', hlrn=HLRN)

fp = int(20 // 0.025) #fixpoint aka t_e

R = np.zeros([2,2,3])
Rstd = np.zeros([2,2,3])

start = int(startps * num2 / maxps)
end = int(endps * num2 / maxps)

print("maxps: ", maxps)
print("num: ", num)
print("num2: ", num2)
print("start: ", start, " end: ", end)

def printStd(a,b):
    print("%.4f +/- %.4f" %(a,b))

N, l1, l2, delL, c1, dc1, c2, dc2, Time2 = stats.CalcSolution(fp, State, TIMESTEPS, start, end, R, Rstd, TRate, TRAJ_COUNTER, num2)


tR = stats.ResTime(l1, N[0][fp], N[1][fp], R[:,:,0])
tRstd = stats.StdResTime(l1, delL, N[0][fp], N[1][fp], R[:,:,0])
tRsimple = 1./math.fabs(l1) * 6.53
U = 90 * 1e-3 * e0
tP = tR * math.exp(-U / (kB*Temp))
print('\ntR:')
printStd(tR, tRstd)
print('\n')

horizontal = np.zeros([4,len(Time)])
horizontal[0] = np.array([R[1,0,0] for i in range(len(Time))])
horizontal[1] = np.array([R[0,1,0] for i in range(len(Time))])
horizontal[2] = np.array([-R[0,0,0]-R[1,0,0] for i in range(len(Time))])
horizontal[3] = np.array([-R[1,1,0]-R[0,1,0] for i in range(len(Time))])

left = State[0,-1000] * horizontal[0,-1000]
right = State[1,-1000] * horizontal[2,-1000] + State[1,-1000] * horizontal[3,-1000]
print('N_T * T_QT \t\tvs.\t N_Q * (T_TQ + T_CQ)')
print(left, '\tvs.\t', right)

#fig, ax = plt.subplots(num='a'+angle+'t'+temperature+'e'+Energy[energy])

plot.TransitionRate(angle, temperature, Energy[energy], TimePrime2, TRate[0, :num2], TRate[2, :num2], TRate[1, :num2], TRate[3, :num2], lblA='T_QT', lblB='T_TQ', lblC='T_CT', lblD='T_CQ',
smflag=0, pltflag=1, ylbl='T / t\u2080\u207B\u00B9', avgflag=1, start=start, end=end, hlrn=HLRN)

#fl = open("ResTime.dat", 'a')
#fl.write("%2d %3d %.3f %.4f %.4f %.4f %.4f\n" %(int(angle), int(temperature), float(Energy[energy]), tR, tRstd, tRsimple/tR, tP))
#fl.close()
#
if HLRN == 1:
    nameA = 'a'+angle+'t'+Temperature[temperature]+'e'+Energy[energy]+'HLRN'
else:
    nameA = 'a'+angle+'t'+Temperature[temperature]+'e'+Energy[energy]

plt.plot(Time2*0.025, N[0], 'b', label='T_an', alpha=0.5)
plt.plot(TimePrime*0.025, StateComp[0], 'b-.', label='T_num')
plt.plot(Time2*0.025, N[1], 'g', label='Q_an', alpha=0.5)
plt.plot(TimePrime*0.025, StateComp[1], 'g--', label='Q_num')
plt.legend()
plt.title("Method 1, Angle " + angle + " deg, Energy " + Energy[energy] + " meV, Temp " + temperature + " K")
plt.axis([0, 60, -.05, 1.05])
plt.tight_layout()
#plt.savefig('/home/becker/lammps/tmpplot/' + nameA + 'Solution.pdf')
plt.show(block=True)
plt.close()

#sys.exit()

# arrays for mean energy value, needed for equilibration analysis
paraAvg = np.zeros([2,30])
normAvg = np.zeros([2,30])
kinAvg  = np.zeros([2,30])
wf = 0
pf = 0
filename = in_folder + 'E' + '.dat'
for n in range(0,30,5):
    #sprint(str(n))
    delta_T = np.zeros([NUM_OF_JOBS*NUM_OF_TRAJ])
    delta_Q = np.zeros([NUM_OF_JOBS*NUM_OF_TRAJ])
    delta_C = np.zeros([NUM_OF_JOBS*NUM_OF_TRAJ])
    E_T = np.zeros([5, NUM_OF_JOBS*NUM_OF_TRAJ])
    E_Q = np.zeros([5, NUM_OF_JOBS*NUM_OF_TRAJ])
    E_C = np.zeros([5, NUM_OF_JOBS*NUM_OF_TRAJ])
    #TODO
    #filename = in_folder + 'deltaE' + str(n) + '.dat'
    t,q,c = nrg.EnergyLoss(n, delta_T, delta_Q, delta_C, NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[4:7,:,:], Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    plot.Histogram(-130,150,141, delta_T[:t], delta_Q[:q], delta_C[:c], subplts=1, lblA='delta T', lblB='delta Q', lblC='delta C', lbl='Total delta E',
    pltflag=pf, sm_flag=1, nu=2, writeflag=wf, flname=filename, nb=n, action=0, scl=1e1)
    #nrg.PlotHist(-120,140,80, delta_T[:t], delta_Q[:q], delta_C[:c], subplts=1, lblA='delta T', lblB='delta Q', lblC='delta C', lbl='Total delta E', pltflag=1, sm_flag=1)

    t,q,c = nrg.kinEnergyDistr(n, E_T[4,:], E_Q[4,:], E_C[4,:], NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    *waste, average, stdev = plot.Histogram(0,1000,100, E_T[4][:t], E_Q[4][:q], E_C[4][:c], lblA='T', lblB='Q', lblC='C', lbl='Kinetic Energy Distr',
    edge='pos', edgeA='pos', edgeB='pos', edgeC='pos', Title='KE, nb= '+str(n),
    subplts=0, sm_flag=1, pltflag=pf, nu=2, writeflag=wf, flname=filename, avg=1, std=1, nb=n, action=5, scl=1e1)
    kinAvg[0,n] = average
    kinAvg[1,n] = stdev

    #filename = in_folder + 'paraE' + str(n) + '.dat'
    t,q,c = nrg.paraEnergyDistr(n, E_T[0,:], E_Q[0,:], E_C[0,:], NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    *waste, average, stdev = plot.Histogram(0,80,80, E_T[0][:t], E_Q[0][:q], E_C[0][:c], lblA='T', lblB='Q', lblC='C', lbl='Parallel Energy Distr',
    edge='pos', edgeA='pos', edgeB='pos', edgeC='pos',
    subplts=1, sm_flag=1, pltflag=pf, nu=1+math.sin(n/30*math.pi/2), writeflag=wf, flname=filename, avg=1, std=1, nb=n, action=1, scl=1e1)
    # with the sine function we increase the smoothing with growing n, as then the number of existing trajectories decreases
    paraAvg[0,n] = average
    paraAvg[1,n] = stdev    # so far in plot.Histogram() we define stdev as (<x^2> - <x>^2)/N. don't know, if it's true
                            # somehow we need to incorporate the lack of trajectories at later bounces
    #filename = in_folder + 'normE' + str(n) + '.dat'
    t,q,c = nrg.normEnergyDistr(n, E_T[1,:], E_Q[1,:], E_C[1,:], NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    *waste, average, stdev = plot.Histogram(0,200,100, E_T[1][:t], E_Q[1][:q], E_C[1][:c], subplts=1, lblA='T', lblB='Q', lblC='C', lbl='Normal Energy Distr',
    edge='pos', edgeA='pos', edgeB='pos', edgeC='pos', avg=1, std=1,
    sm_flag=1, pltflag=pf, nu=2, writeflag=wf, flname=filename, nb=n, action=2, scl=1e1)
    normAvg[0,n] = average
    normAvg[1,n] = stdev    # so far in plot.Histogram() we define stdev as (<x^2> - <x>^2)/N. don't know, if it's true
                            # somehow we need to incorporate the lack of trajectories at later bounces

    #filename = in_folder + 'potE' + str(n) + '.dat'
    t,q,c = nrg.potEnergyDistr(n, E_T[2,:], E_Q[2,:], E_C[2,:], NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    plot.Histogram(-140,0,200, E_T[2][:t], E_Q[2][:q], E_C[2][:c], subplts=1, lblA='T', lblB='Q', lblC='C', lbl='Potential Energy Distr', sm_flag=1, pltflag=pf, nu=1, writeflag=wf, flname=filename,
    nb=n, action=3, scl=1e1)

    #filename = in_folder + 'totE' + str(n) + '.dat'
    t,q,c = nrg.totEnergyDistr(n, E_T[3,:], E_Q[3,:], E_C[3,:], NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    plot.Histogram(-180,180,145, E_T[3][:t], E_Q[3][:q], E_C[3][:c], subplts=1, lblA='T', lblB='Q', lblC='C', lbl='Total Energy Distr', edge='none', edgeA='neg', edgeB='pos', edgeC='pos',
    sm_flag=1, pltflag=pf, nu=2, writeflag=wf, flname=filename, nb=n, action=4, scl=1e1)


#evaluate how many bounces typically appear in a trajectory with fixed timesteps at Nt = 120000
sequence = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
bouncehist = np.histogram(bounceavg, density=False, bins=sequence)

#look at how the number of trajectories diminishes as nb grows
bouncecumulative = np.cumsum(bouncehist[0])

#calc StdDev for energy development over time
for n in range(0,1):
    paraAvg[1,n] = stats.stdDev(Bounce[9,n,:], paraAvg[0,n], (NUM_OF_JOBS*NUM_OF_TRAJ)-bouncecumulative[n])
    normAvg[1,n] = stats.stdDev(Bounce[10,n,:], normAvg[0,n], (NUM_OF_JOBS*NUM_OF_TRAJ)-bouncecumulative[n])

#calc average bounce time + stddev
TimeArr = np.zeros([3,80])
TimeArr = bounce.BncTime(Bounce, TOTAL, nbnc, TIMESTEPS, TimeArr)
timebetween = bounce.TimeBetwBounces(TimeArr[:,:nbnc]) * 2.5e-2
print(timebetween)
plt.errorbar(TimeArr[0][:nbnc]*2.5e-2, np.arange(0,nbnc,1), xerr=TimeArr[1][:nbnc]*2.5e-2)
plt.xlabel("Time / ps")
plt.ylabel("Bounces")
plt.show()
