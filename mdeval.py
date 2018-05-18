#mdeval.py
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import plot
import bounce
import nrg
import smooth
import stats

Angle  = {'0':'0.00', '30':'0.52', '45':'0.80', '60':'1.05'}
Temperature={'80':'80', '190':'190', '300':'300'}
try:
    angle = sys.argv[1]
    temperature = sys.argv[2]
except:
    angle = '30'
    temperature = '300'
    print('Opening default: a0.52t300e16.027')
try:
    energy = sys.argv[3]
    Energy = {'12':'12.95', '16':'16.027', '21':'21.71', '36':'36.55', '70':'70.0'}
except:
    energy = angle
    Energy = {'0':'12.95', '30':'16.027', '45':'21.71', '60':'36.55', '70':'70.0'}
jobs = (7227490, 7227494, 7227495, 7227496, 7227497)
if energy != '70':
    if int(temperature) == 80:
        jobs = (7094790, 7107184, 7107185, 7107186, 7107187, 7107195, 7131096, 7131097, 7131115, 7131116)
        #jobs = (7279669, 7279673, 7279675, 7279676, 7279677)
        te = 11 // 0.025 #T=80 K
    elif int(temperature) == 190:
        jobs = (7237382, 7237376, 7237379, 7237380, 7237381)
        te = 15 // 0.025 #T=190 K
    elif int(temperature) == 300:
        jobs = (7094785, 7107197, 7107198, 7107199, 7107200, 7107201, 7131101, 7131102, 7131113, 7131114)
        te = 15 // 0.025 #T=300 K
name = 'a'+Angle[angle]+'t'+Temperature[temperature]+'e'+Energy[energy]
''' Ar-Pt analysis
jobs = (7331694, 7331705, 7331706, 7331713, 7331723)
name = 'Pta'+Angle[angle]+'t'+Temperature[temperature]+'e'+Energy[energy]
'''
in_folder = '/home/becker/lammps/111/' + name + '/'

TIMESTEPS = 120000 // 100
DT = 0.00025
NUM_OF_TRAJ = 200
NUM_OF_JOBS = len(jobs)
#SCALING = 1.0/(NUM_OF_JOBS*NUM_OF_TRAJ)
SCALING = 1
TOTAL = np.zeros([TIMESTEPS])
Time = np.arange(0,TIMESTEPS,1)
State = np.zeros([3,TIMESTEPS]) # T, Q, C
Trans = np.zeros([6,TIMESTEPS]) # QT, CT, TQ, CQ, TC, QC
TRate = np.zeros([6,TIMESTEPS]) # QT, CT, TQ, CQ, TC, QC
Refl = np.zeros([TIMESTEPS])
Bounce = np.zeros([13,80,NUM_OF_JOBS*NUM_OF_TRAJ])
# nb, tm, s1, t1, E_xy.i, E_z.i, V.i, s2, t2, E_xy.f, E_z.f, V.f, refl_flag

TOTAL = NUM_OF_JOBS*NUM_OF_TRAJ
bounceavg = np.zeros([TOTAL])

print('Evaluating Data for', name)
for ctr, jb in enumerate(jobs):
    for d in range(1, NUM_OF_TRAJ+1):
        flag = 0
        # allocate Trajectory array
        # Ez, Exy, V, vz
        Traj = np.zeros([4,TIMESTEPS])

        # open file for potential energy;
        # keep in mind that relevant information is only stored in every tenth line
        name_pe = in_folder + str(jb) + '/' + str(d) + ".lammpstrj"
        fl_pe = open(name_pe, 'r')
        j = 1
        for i, line in enumerate(fl_pe):
            assert(j <= 10)
            if j == 10:
                idx = i//10
                epot = float(line.split()[0])
                Traj[2, idx-1] = epot*2.0
                j = 1
            else:
                j += 1
        #end for (potential file)
        # open file with trajectory information;
        # here except for the header every line contains information for the
        # corresponding timestep
        name_kin = in_folder + str(jb) + '/' + str(d) + ".dat"
        fl_kin = open(name_kin, 'r')
        for a,line in enumerate(fl_kin):
            b = a - 1
            if a == 0:
                continue

            t, x, y, z, vx, vy, Traj[3,b], en = line.split()
            Traj[0,b] = bounce.normalKE(float(Traj[3,b]))
            Traj[1,b] = bounce.parallelKE(float(vx), float(vy))
            assert(Traj[0,b] >= 0)
            assert(Traj[1,b] >= 0)
        #end for (kinetic file)
        jj = ctr*200 + d - 1
        nb = bounce.Countbounces(Traj[3,:], Traj[1,:], Traj[0,:], Traj[2,:], Bounce[:,:,jj],jj)
        if nb == -9999:
            TOTAL -= 1
            continue

        for i in range(0,nb):
            try:
                flag = bounce.Population(i, Bounce[:,i-1,jj], Bounce[:,i,jj], Bounce[:,i+1,jj], State, Trans, Refl, 1., jj)
            except:
                flag = bounce.Population(i, Bounce[:,i,jj], Bounce[:,i,jj], Bounce[:,i+1,jj], State, Trans, Refl, 1., jj)
            if flag == -6000 or flag == -5000:
                TOTAL -= 1
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
        fl_kin.close()
        fl_pe.close()

    #end for (trajectory)
#end for (job list)

State[0,:] /= TOTAL
State[1,:] /= TOTAL
State[2,:] /= TOTAL

Trans[0,:] /= TOTAL
Trans[1,:] /= TOTAL
Trans[2,:] /= TOTAL
Trans[3,:] /= TOTAL
Trans[4,:] /= TOTAL
Trans[5,:] /= TOTAL
Refl[:] /= TOTAL

StateNew = np.zeros([3, TIMESTEPS])
StateSum = np.zeros([3, TIMESTEPS])
TransNew = np.zeros([6, TIMESTEPS])
TRateE = np.zeros([6, TIMESTEPS])
TRateA = np.zeros([6, TIMESTEPS])
                #4       2              0           1
StateSum[0] = Trans[4] + Trans[2] - Trans[0] - Trans[1]
                #5            0           2           3
StateSum[1] = Trans[5] + Trans[0] - Trans[2] - Trans[3]
#C              CT          CQ          TC      QC
StateSum[2] = 1 + Trans[1] + Trans[3] - Trans[4] - Trans[5]

StickCoeff = stats.InitSticking(Refl, float(Angle[angle]),float(Temperature[temperature]), float(Energy[energy]), TOTAL)

plot.Populations(angle, temperature, Energy[energy], Time, State[0], State[1], State[2], smflag=0, pltflag=0)
Trans[0], Trans[2], Trans[1], Trans[3], Trans[4], Trans[5] = plot.TransitionPopulations(angle, temperature, Energy[energy], Time, Trans[0], Trans[2], Trans[1], Trans[3], np.zeros([TIMESTEPS]),np.zeros([TIMESTEPS]), smflag=0, pltflag=0)

'''
print("Plot Summed Populations")
fig, ax = plt.subplots()
ax.plot(Time*2.5e-2, StateSum[0], '-', label = 'Trapped')
ax.plot(Time*2.5e-2, State[0], '-', label = 'TrappedNum')
ax.plot(Time*2.5e-2, StateSum[1], '-', label = 'Quasi')
ax.plot(Time*2.5e-2, State[1], '-', label = 'QuasiNum')
ax.plot(Time*2.5e-2, StateSum[2], '-', label = 'Cont')
ax.plot(Time*2.5e-2, State[2], '-', label = 'ContNum')
ax.plot(Time*2.5e-2, StateSum[0] + StateSum[1] + StateSum[2], '-', label = 'Norm check')
ax.plot(Time*2.5e-2, State[0] + State[1] + State[2], '-', label = 'Norm check Num')
legend = ax.legend()
plt.xlabel("time / ps")
plt.ylabel("")
plt.title("Comparison of numerical and summed populations")
plt.show(block=True)
plt.close(fig)
'''

'''State = np.zeros([3,TIMESTEPS]) # T, Q, C
Trans = np.zeros([4,TIMESTEPS]) # QT, CT, TQ, CQ
TRate = np.zeros([4, TIMESTEPS]) # QT, CT, TQ, CQ, TC, QC'''

#TRate[0] = bounce.TransitionRate(1, State[0], Trans[0], TRate[0]) # QT
#TRate[1] = bounce.TransitionRate(1, State[0], Trans[1], TRate[1]) # CT
#TRate[2] = bounce.TransitionRate(1, State[1], Trans[2], TRate[2]) # TQ
#TRate[3] = bounce.TransitionRate(1, State[1], Trans[3], TRate[3]) # CQ
#TRate[4] = bounce.TransitionRate(1, State[2], Trans[4], TRate[4]) # TC
#TRate[5] = bounce.TransitionRate(1, State[2], Trans[5], TRate[5]) # QC

print("Calculate Transition Rates")
d = 0
for i in range(0,6):
    TRate[i] = bounce.DifferentiateT(30, State[d], Trans[i], TRate[i])
    d += i % 2

d = 0
for i in range(0,6):
    TRateE[i] = bounce.TransitionRateE(Trans[i], State[d], te, TRateE[i], TIMESTEPS)
    d += i % 2


d = 0
for i in range(0,6):
    TRateA[i] = bounce.AlternTransition(Trans[i], State[d], int(10//0.025), TRateA[i], TIMESTEPS)
    d += i % 2


left = State[0,-400] * TRate[0,-400]
right = State[1,-400] * TRate[2,-400] + State[1,-400] * TRate[3,-400]
print('N_T * T_QT \t\tvs.\t N_Q * (T_TQ + T_CQ)')
print(left, '\tvs.\t', right)
plot.TransitionRate(angle, temperature, Energy[energy], Time, TRate[0], TRate[2], TRate[1], TRate[3], lblA='T_QT', lblB='T_TQ', lblC='T_CT', lblD='T_CQ', smflag=0, pltflag=0, ylbl='T / t\u2080\u207B\u00B9')
plot.TransitionRate(angle, temperature, Energy[energy], Time, TRateE[0], TRateE[2], TRateE[1], TRateE[3], lblA='T_QT', lblB='T_TQ', lblC='T_CT', lblD='T_CQ', smflag=0, pltflag=0)
plot.TransitionRate(angle, temperature, Energy[energy], Time, TRateA[0], TRateA[2], TRateA[1], TRateA[3], lblA='T_QT', lblB='T_TQ', lblC='T_CT', lblD='T_CQ', smflag=0, pltflag=0, ylbl='Alternative Calculation')

StateNew[0] = bounce.IntegratePopulationT(StateNew[0], State, TRate, 1)
StateNew[1] = bounce.IntegratePopulationQ(StateNew[1], State, TRate, 1)
StateNew[2] = bounce.IntegratePopulationC(StateNew[2], State, TRate, 1)
plot.Populations(angle, temperature, Energy[energy], Time, StateNew[0], StateNew[1], StateNew[2], smflag=0, pltflag=0)
plot.Populations(angle, temperature, Energy[energy], Time, StateSum[0], StateSum[1], StateSum[2], smflag=0, pltflag=0)
#sys.exit()
R = np.zeros([2,2,3])
Rstd = np.zeros([2,2,3])
start = int(15//0.025)
end = int(30//0.025)
stats.CalcR(start, end, R[:,:,0], Rstd[:,:,0], TRate)

# arrays for mean energy value, needed for equilibration analysis
paraAvg = np.zeros([2,30])
normAvg = np.zeros([2,30])
wf = 0
filename = in_folder + 'E' + '.dat'
for n in range(0,30):

    #sprint(str(n))
    delta_T = np.zeros([NUM_OF_JOBS*NUM_OF_TRAJ])
    delta_Q = np.zeros([NUM_OF_JOBS*NUM_OF_TRAJ])
    delta_C = np.zeros([NUM_OF_JOBS*NUM_OF_TRAJ])
    E_T = np.zeros([4, NUM_OF_JOBS*NUM_OF_TRAJ])
    E_Q = np.zeros([4, NUM_OF_JOBS*NUM_OF_TRAJ])
    E_C = np.zeros([4, NUM_OF_JOBS*NUM_OF_TRAJ])
    #TODO
    #filename = in_folder + 'deltaE' + str(n) + '.dat'
    t,q,c = nrg.EnergyLoss(n, delta_T, delta_Q, delta_C, NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[4:7,:,:], Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    plot.Histogram(-130,150,141, delta_T[:t], delta_Q[:q], delta_C[:c], subplts=1, lblA='delta T', lblB='delta Q', lblC='delta C', lbl='Total delta E',
    pltflag=0, sm_flag=1, nu=2, writeflag=wf, flname=filename, nb=n, action=0, scl=1e1)
    #nrg.PlotHist(-120,140,80, delta_T[:t], delta_Q[:q], delta_C[:c], subplts=1, lblA='delta T', lblB='delta Q', lblC='delta C', lbl='Total delta E', pltflag=1, sm_flag=1)

    #filename = in_folder + 'paraE' + str(n) + '.dat'
    t,q,c = nrg.paraEnergyDistr(n, E_T[0,:], E_Q[0,:], E_C[0,:], NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    *waste, average, stdev = plot.Histogram(0,80,80, E_T[0][:t], E_Q[0][:q], E_C[0][:c], lblA='T', lblB='Q', lblC='C', lbl='Parallel Energy Distr',
    edge='pos', edgeA='pos', edgeB='pos', edgeC='pos',
    subplts=1, sm_flag=1, pltflag=0, nu=1+math.sin(n/30*math.pi/2), writeflag=wf, flname=filename, avg=1, std=1, nb=n, action=1, scl=1e1)
    # with the sine function we increase the smoothing with growing n, as then the number of existing trajectories decreases
    paraAvg[0,n] = average
    paraAvg[1,n] = stdev    # so far in plot.Histogram() we define stdev as (<x^2> - <x>^2)/N. don't know, if it's true
                            # somehow we need to incorporate the lack of trajectories at later bounces
    #filename = in_folder + 'normE' + str(n) + '.dat'
    t,q,c = nrg.normEnergyDistr(n, E_T[1,:], E_Q[1,:], E_C[1,:], NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    *waste, average, stdev = plot.Histogram(0,200,100, E_T[1][:t], E_Q[1][:q], E_C[1][:c], subplts=1, lblA='T', lblB='Q', lblC='C', lbl='Normal Energy Distr',
    edge='pos', edgeA='pos', edgeB='pos', edgeC='pos', avg=1, std=1,
    sm_flag=1, pltflag=0, nu=2, writeflag=1, flname=filename, nb=n, action=2, scl=1e1)
    normAvg[0,n] = average
    normAvg[1,n] = stdev    # so far in plot.Histogram() we define stdev as (<x^2> - <x>^2)/N. don't know, if it's true
                            # somehow we need to incorporate the lack of trajectories at later bounces

    #filename = in_folder + 'potE' + str(n) + '.dat'
    t,q,c = nrg.potEnergyDistr(n, E_T[2,:], E_Q[2,:], E_C[2,:], NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    plot.Histogram(-140,0,200, E_T[2][:t], E_Q[2][:q], E_C[2][:c], subplts=1, lblA='T', lblB='Q', lblC='C', lbl='Potential Energy Distr', sm_flag=1, pltflag=0, nu=1, writeflag=wf, flname=filename,
    nb=n, action=3, scl=1e1)

    #filename = in_folder + 'totE' + str(n) + '.dat'
    t,q,c = nrg.totEnergyDistr(n, E_T[3,:], E_Q[3,:], E_C[3,:], NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    plot.Histogram(-180,180,145, E_T[3][:t], E_Q[3][:q], E_C[3][:c], subplts=1, lblA='T', lblB='Q', lblC='C', lbl='Total Energy Distr', edge='none', edgeA='neg', edgeB='pos', edgeC='pos',
    sm_flag=1, pltflag=0, nu=2, writeflag=wf, flname=filename, nb=n, action=4, scl=1e1)




#evaluate how many bounces typically appear in a trajectory with fixed timesteps at Nt = 120000
sequence = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
bouncehist = np.histogram(bounceavg, density=False, bins=sequence)

#look at how the number of trajectories diminishes as nb grows
bouncecumulative = np.cumsum(bouncehist[0])
#plt.plot(np.arange(0,39,1), (NUM_OF_JOBS*NUM_OF_TRAJ)-bouncecumulative)
#plt.xlabel("Bounces")
#plt.ylabel("Trajectories")
#plt.show()
'''
'''
#calc StdDev for energy development over time
for n in range(0,30):
    paraAvg[1,n] = stats.stdDev(Bounce[9,n,:], paraAvg[0,n], (NUM_OF_JOBS*NUM_OF_TRAJ)-bouncecumulative[n])
    normAvg[1,n] = stats.stdDev(Bounce[10,n,:], normAvg[0,n], (NUM_OF_JOBS*NUM_OF_TRAJ)-bouncecumulative[n])

'''
#calc average bounce time + stddev
TimeArr = np.zeros([3,80])
TimeArr = bounce.BncTime(Bounce, TOTAL, 30, TIMESTEPS, TimeArr)

plt.errorbar(TimeArr[0][:30]*2.5e-2, np.arange(0,30,1), xerr=TimeArr[1][:30]*2.5e-2)
plt.xlabel("Time / ps")
plt.ylabel("Bounces")
plt.show()


filename = in_folder + 'Bnc.dat'
fbnc = open(filename, 'w')
fbnc.write("# nb tm delta_tm num of trajectories\n")

for i in range(0,30):
    fbnc.write("%e %e %e %e\n" %(i, TimeArr[0,i]*2.5e-2, TimeArr[1,i]*2.5e-2, TimeArr[2,i]))
fbnc.close()


print('Parallel Component:')
for i in range(0,30):
    print(str(paraAvg[0,i]) + ' +/- ' + str(paraAvg[1,i]))

print('\nNormal Component:')
for i in range(0,30):
    print(str(normAvg[0,i]) + ' +/- ' + str(normAvg[1,i]))


plt.plot(np.arange(0,30,1), paraAvg[0])
plt.show()
plt.plot(np.arange(0,30,1), normAvg[0])
plt.show()
'''

filename = in_folder + 'paraEofBnc.dat'
BncEpara = open(filename, 'w')
BncEpara.write('# nb E deltaE\n')
for i in range(0,30):
    BncEpara.write("%e %e %e\n" %(i, paraAvg[0,i], paraAvg[1,i]))
BncEpara.close()

filename = in_folder + 'normEofBnc.dat'
BncEnorm = open(filename, 'w')
BncEnorm.write('# nb E deltaE\n')
for i in range(0,30):
    BncEnorm.write("%e %e %e\n" %(i, normAvg[0,i], normAvg[1,i]))
BncEnorm.close()
