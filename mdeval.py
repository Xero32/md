#mdeval.py
import numpy as np
import matplotlib.pyplot as plt
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
    elif int(temperature) == 190:
        jobs = (7237382, 7237376, 7237379, 7237380, 7237381)
    elif int(temperature) == 300:
        jobs = (7094785, 7107197, 7107198, 7107199, 7107200, 7107201, 7131101, 7131102, 7131113, 7131114)

name = 'a'+Angle[angle]+'t'+Temperature[temperature]+'e'+Energy[energy]
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
        j = ctr*200 + d - 1
        nb = bounce.Countbounces(Traj[3,:], Traj[1,:], Traj[0,:], Traj[2,:], Bounce[:,:,j])
        if nb == -9999:
            TOTAL -= 1
            continue

        for i in range(0,nb):
            flag = bounce.Population(i, Bounce[:,i,j], Bounce[:,i+1,j], State, Trans, Refl, 1., j)
            if flag == -6000 or flag == -5000:
                TOTAL -= 1
                print(str(flag))
                break
            elif flag == -3000:
                break
            if flag == 1:
                Bounce[12,i,j] = flag
                break

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
                #4       2              0           1
StateSum[0] = Trans[4] + Trans[2] - Trans[0] - Trans[1]
                #5            0           2           3
StateSum[1] = Trans[5] + Trans[0] - Trans[2] - Trans[3]
#C              CT          CQ          TC      QC
StateSum[2] = 1 + Trans[1] + Trans[3] - Trans[4] - Trans[5]

StickCoeff = stats.InitSticking(Refl, float(Angle[angle]),float(Temperature[temperature]), float(Energy[energy]), TOTAL)

plot.Populations(angle, temperature, Energy[energy], Time, State[0], State[1], State[2], smflag=1, pltflag=1)
plot.TransitionPopulations(angle, temperature, Energy[energy], Time, Trans[0], Trans[2], Trans[1], Trans[3], Trans[4], Trans[5], smflag=1, pltflag=1)

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


'''State = np.zeros([3,TIMESTEPS]) # T, Q, C
Trans = np.zeros([4,TIMESTEPS]) # QT, CT, TQ, CQ
TRate = np.zeros([4, TIMESTEPS]) # QT, CT, TQ, CQ, TC, QC'''

TRate[0] = bounce.TransitionRate(1, StateSum[0], Trans[0], TRate[0]) # QT
TRate[2] = bounce.TransitionRate(1, StateSum[1], Trans[2], TRate[2]) # TQ
TRate[3] = bounce.TransitionRate(1, StateSum[1], Trans[3], TRate[3]) # CQ
TRate[1] = bounce.TransitionRate(1, StateSum[0], Trans[1], TRate[1]) # CT
TRate[4] = bounce.TransitionRate(1, StateSum[2], Trans[4], TRate[4]) # TC
TRate[5] = bounce.TransitionRate(1, StateSum[2], Trans[5], TRate[5]) # QC
left = State[0,-400] * TRate[0,-400]
right = State[1,-400] * TRate[2,-400] + State[1,-400] * TRate[3,-400]
print('N_T * T_QT \t\tvs.\t N_Q * (T_TQ + T_CQ)')
print(left, '\tvs.\t', right)
plot.TransitionRate(angle, temperature, Energy[energy], Time, TRate[0], TRate[2], TRate[1], TRate[3], TRate[4], TRate[5], lblA='T_QT', lblB='T_TQ', lblC='T_CT', lblD='T_CQ', smflag=1, pltflag=1)

StateNew[0] = bounce.IntegratePopulationT(StateNew[0], StateSum, TRate, 1)
StateNew[1] = bounce.IntegratePopulationQ(StateNew[1], StateSum, TRate, 1)
StateNew[2] = bounce.IntegratePopulationC(StateNew[2], StateSum, TRate, 1)
plot.Populations(angle, temperature, Energy[energy], Time, StateNew[0], StateNew[1], StateNew[2], smflag=0, pltflag=0)
print("Plot Summed Populations")
fig, ax = plt.subplots()
ax.plot(Time*2.5e-2, StateNew[0], '-', label = 'Trapped')
ax.plot(Time*2.5e-2, State[0], '-', label = 'TrappedNum')
ax.plot(Time*2.5e-2, StateNew[1], '-', label = 'Quasi')
ax.plot(Time*2.5e-2, State[1], '-', label = 'QuasiNum')
ax.plot(Time*2.5e-2, StateNew[2], '-', label = 'Cont')
ax.plot(Time*2.5e-2, State[2], '-', label = 'ContNum')
ax.plot(Time*2.5e-2, StateNew[0] + StateNew[1] + StateNew[2], '-', label = 'Norm check')
ax.plot(Time*2.5e-2, State[0] + State[1] + State[2], '-', label = 'Norm check Num')
legend = ax.legend()
plt.xlabel("time / ps")
plt.ylabel("")
plt.title("Comparison of numerical and integrated populations")
plt.show(block=True)
plt.close(fig)

'''
for n in range(0,30):
    #sprint(str(n))
    delta_T = np.zeros([NUM_OF_JOBS*NUM_OF_TRAJ])
    delta_Q = np.zeros([NUM_OF_JOBS*NUM_OF_TRAJ])
    delta_C = np.zeros([NUM_OF_JOBS*NUM_OF_TRAJ])
    E_T = np.zeros([4, NUM_OF_JOBS*NUM_OF_TRAJ])
    E_Q = np.zeros([4, NUM_OF_JOBS*NUM_OF_TRAJ])
    E_C = np.zeros([4, NUM_OF_JOBS*NUM_OF_TRAJ])
    #TODO

    filename = in_folder + 'deltaE' + str(n) + '.dat'
    t,q,c = nrg.EnergyLoss(n, delta_T, delta_Q, delta_C, NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[4:7,:,:], Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    plot.Histogram(-120,140,80, delta_T[:t], delta_Q[:q], delta_C[:c], subplts=1, lblA='delta T', lblB='delta Q', lblC='delta C', lbl='Total delta E', pltflag=0, sm_flag=1, nu=1, writeflag=0, flname=filename)
    #nrg.PlotHist(-120,140,80, delta_T[:t], delta_Q[:q], delta_C[:c], subplts=1, lblA='delta T', lblB='delta Q', lblC='delta C', lbl='Total delta E', pltflag=1, sm_flag=1)

    filename = in_folder + 'paraE' + str(n) + '.dat'
    t,q,c = nrg.paraEnergyDistr(n, E_T[0,:], E_Q[0,:], E_C[0,:], NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    plot.Histogram(-5,60,30, E_T[0][:t], E_Q[0][:q], E_C[0][:c], subplts=1, lblA='T', lblB='Q', lblC='C', lbl='Parallel Energy Distr', sm_flag=1, pltflag=0, edge='pos', edgeA='pos', edgeB='pos', edgeC='pos', nu=1, writeflag=0, flname=filename)

    filename = in_folder + 'normE' + str(n) + '.dat'
    t,q,c = nrg.normEnergyDistr(n, E_T[1,:], E_Q[1,:], E_C[1,:], NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    plot.Histogram(-5,200,41, E_T[1][:t], E_Q[1][:q], E_C[1][:c], subplts=1, lblA='T', lblB='Q', lblC='C', lbl='Normal Energy Distr', sm_flag=1, pltflag=0, edge='pos', edgeA='pos', edgeB='pos', edgeC='pos', nu=1, writeflag=0, flname=filename)

    filename = in_folder + 'potE' + str(n) + '.dat'
    t,q,c = nrg.potEnergyDistr(n, E_T[2,:], E_Q[2,:], E_C[2,:], NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    plot.Histogram(-140,0,100, E_T[2][:t], E_Q[2][:q], E_C[2][:c], subplts=1, lblA='T', lblB='Q', lblC='C', lbl='Potential Energy Distr', sm_flag=1, pltflag=0, nu=1, writeflag=1, flname=filename)

    filename = in_folder + 'totE' + str(n) + '.dat'
    t,q,c = nrg.totEnergyDistr(n, E_T[3,:], E_Q[3,:], E_C[3,:], NUM_OF_JOBS*NUM_OF_TRAJ, Bounce[9:12,:,:], Bounce[7,:,:], Bounce[12,:,:])
    plot.Histogram(-180,180,100, E_T[3][:t], E_Q[3][:q], E_C[3][:c], subplts=1, lblA='T', lblB='Q', lblC='C', lbl='Total Energy Distr', sm_flag=1, pltflag=0, edge='none', edgeA='neg', edgeB='pos', edgeC='pos', nu=1, writeflag=0, flname=filename)
'''
