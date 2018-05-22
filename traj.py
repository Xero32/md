#traj.py
import numpy as np
import matplotlib.pyplot as plt
import sys
import plot
import bounce
import nrg
import smooth
import stats
TIMESTEPS = 525000 // 500
DT = 0.00006
NUM_OF_TRAJ = 200
#NUM_OF_JOBS = len(jobs)
NUM_OF_JOBS = 10
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
Traj = np.zeros([5,TIMESTEPS])
name_pe = './111/a0.52t80e16.027/LONG/7275348/58.lammpstrj'
fl_pe = open(name_pe, 'r')
j = 1
ctr = 0
d = 11
for i, line in enumerate(fl_pe):
    assert(j <= 50)
    if j == 50:
        idx = i//50
        epot = float(line.split()[0])
        Traj[2, idx-1] = epot*2.0
        j = 1
    else:
        j += 1
#end for (potential file)
# open file with trajectory information;
# here except for the header every line contains information for the
# corresponding timestep
7094785
name_kin = './111/a0.52t80e16.027/LONG/7275348/58
.dat'
fl_kin = open(name_kin, 'r')
for a,line in enumerate(fl_kin):
    b = a - 1
    if a == 0:
        continue

    t, x, y, Traj[4,b], vx, vy, Traj[3,b], en = line.split()
    Traj[0,b] = bounce.normalKE(float(Traj[3,b]))
    Traj[1,b] = bounce.parallelKE(float(vx), float(vy))
    assert(Traj[0,b] >= 0)
    assert(Traj[1,b] >= 0)
#end for (kinetic file)
jj = ctr*200 + d - 1
nb = bounce.Countbounces(Traj[3,:], Traj[1,:], Traj[0,:], Traj[2,:], Bounce[:,:,jj],jj)


for i in range(0,nb):
    try:
        flag = bounce.Population(i, Bounce[:,i-1,jj], Bounce[:,i,jj], Bounce[:,i+1,jj], State, Trans, Refl, 1., jj)
    except:
        flag = bounce.Population(i, Bounce[:,i,jj], Bounce[:,i,jj], Bounce[:,i+1,jj], State, Trans, Refl, 1., jj)
    if flag == -6000 or flag == -5000:
        TOTAL -= 1
        print(str(flag))
        break
    elif flag == -3000:
        break
    if flag == 1:
        Bounce[12,i,jj] = flag
        break

fl_kin.close()
fl_pe.close()

fl = open('trajtestdt.dat','w')
fl.write('# Ez Exy V vz z state\n')
for i in range(0,TIMESTEPS):
    fl.write(str(i) + ' ')
    fl.write(str(Traj[0,i] * 1e3) + ' ')
    fl.write(str(Traj[1,i] * 1e3) + ' ')
    fl.write(str(Traj[2,i] * 1e3) + ' ')
    fl.write(str(Traj[3,i]) + ' ')
    fl.write(str(Traj[4,i]) + ' ')
    fl.write(str(bounce.State(Traj[0,i], Traj[1,i], Traj[2,i])) + '\n')

fl.close()
