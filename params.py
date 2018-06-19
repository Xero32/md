# input.py
import numpy as np
import bounce
import argparse

def InpParams(parser):
    parser.add_argument("a", help='Angle', default='30', nargs='?')
    parser.add_argument("t", help='Temperature', default='300', nargs='?')
    parser.add_argument("--nrg", help='Incident Energy')
    parser.add_argument("--start", help="Sets start time to calculate average transition rate", default=20)
    parser.add_argument("--end", help="Sets end time to calculate average transition rate", default=30)
    parser.add_argument("--hlrn", help="Use trajectories from HLRN runs", action='store_true')
    args = parser.parse_args()
    if args.start:
        startps = float(args.start)
        print(startps)
    if args.end:
        endps = float(args.end)
        print(endps)
    if args.a:
        angle = args.a
    else:
        angle = '30'
    if args.t:
        temperature = args.t
    else:
        temperature = '300'
    if args.nrg:
        energy = args.nrg
    else:
        energy = args.a
    if args.hlrn:
        HLRN = 1
    else:
        HLRN = 0

    return startps, endps, angle, temperature, energy, HLRN

def Parameters(temperature, HLRN):
    if int(temperature) == 80:
        if HLRN:
            TIMESTEPS_RZ = 120000 // 100
            TIMESTEPS_HLRN = 160000 // 100
            NUM_OF_TRAJ_RZ = 200
            NUM_OF_TRAJ = 1000
            nbnc = 40
            jobs = (1931936, 1931937)
        else:
            TIMESTEPS_RZ = 120000 // 100
            TIMESTEPS_HLRN = 120000 // 100
            NUM_OF_TRAJ_RZ = 200
            NUM_OF_TRAJ = 200
            nbnc = 30
            jobs = (7094790, 7107184, 7107185, 7107186, 7107187, 7107195, 7131096, 7131097, 7131115, 7131116)
        #jobs = (7279669, 7279673, 7279675, 7279676, 7279677)
        te = 11 // 0.025 #T=80 K
    elif int(temperature) == 120:
        if HLRN:
            TIMESTEPS_RZ = 120000 // 100
            TIMESTEPS_HLRN = 160000 // 100
            NUM_OF_TRAJ_RZ = 200
            NUM_OF_TRAJ = 1000
            nbnc = 40
            jobs = (1931939, 1935306)
        else:
            TIMESTEPS_RZ = 120000 // 100
            TIMESTEPS_HLRN = 120000 // 100
            NUM_OF_TRAJ_RZ = 200
            NUM_OF_TRAJ = 200
            nbnc = 30
            jobs = ()
    elif int(temperature) == 160:
        if HLRN:
            TIMESTEPS_RZ = 120000 // 100
            TIMESTEPS_HLRN = 160000 // 100
            NUM_OF_TRAJ_RZ = 200
            NUM_OF_TRAJ = 1000
            nbnc = 40
            jobs = (1931943, 1935310)

        else:
            TIMESTEPS_RZ = 120000 // 100
            TIMESTEPS_HLRN = 120000 // 100
            NUM_OF_TRAJ_RZ = 200
            NUM_OF_TRAJ = 200
            nbnc = 30
            jobs = ()
    elif int(temperature) == 190:
        if HLRN:
            TIMESTEPS_RZ = 120000 // 100
            TIMESTEPS_HLRN = 160000 // 100
            NUM_OF_TRAJ_RZ = 200
            NUM_OF_TRAJ = 1000
            nbnc = 40
            jobs = (1935315, 1935316)
        else:
            TIMESTEPS_RZ = 120000 // 100
            TIMESTEPS_HLRN = 120000 // 100
            NUM_OF_TRAJ_RZ = 200
            NUM_OF_TRAJ = 200
            nbnc = 30
            jobs = (7237382, 7237376, 7237379, 7237380, 7237381, 8199182, 8199181, 8199180, 8199179, 8199176)
    elif int(temperature) == 240:
        if HLRN:
            TIMESTEPS_RZ = 120000 // 100
            TIMESTEPS_HLRN = 160000 // 100
            NUM_OF_TRAJ_RZ = 200
            NUM_OF_TRAJ = 1000
            nbnc = 40
            jobs = (1931944, 1931945)
        else:
            TIMESTEPS_RZ = 120000 // 100
            TIMESTEPS_HLRN = 120000 // 100
            NUM_OF_TRAJ_RZ = 200
            NUM_OF_TRAJ = 200
            nbnc = 30
            jobs = ()
    elif int(temperature) == 300:
        TIMESTEPS_RZ = 120000 // 100
        NUM_OF_TRAJ_RZ = 200
        BOTH = 0
        if HLRN == 0:
            TIMESTEPS_HLRN = 120000 // 100
            NUM_OF_TRAJ = 200
            nbnc = 30
            jobs = (7094785, 7107197, 7107198, 7107199, 7107200, 7107201, 7131101, 7131102, 7131113, 7131114)
        else:
            TIMESTEPS_HLRN = 200000 // 100
            NUM_OF_TRAJ = 500
            nbnc = 50
            if BOTH == 1:
                jobs = (7094785, 7107197, 7107198, 7107199, 7107200, 7107201, 7131101, 7131102, 7131113, 7131114, 1922961, 1922962, 1922963, 1922964) # rz + hlrn jobs
            else:
                jobs = (1922961, 1922962, 1922963, 1922964) # hlrn jobs

    return TIMESTEPS_RZ, TIMESTEPS_HLRN, NUM_OF_TRAJ_RZ, NUM_OF_TRAJ, nbnc, jobs

def Openfile(d, in_folder, name, jb, ending, TIMESTEPS_RZ, TIMESTEPS_HLRN, NUM_OF_TRAJ_RZ):
    try:
        try:
            name_kin = in_folder + str(jb) + '.bbatch.hsn.hlrn.de/' + str(d) + ending
            fl_kin = open(name_kin, 'r')
        except:
            name_kin = '/home/becker/lammps/111/HLRN/' + name + '/' + str(jb) + '.bbatch.hsn.hlrn.de/' + str(d) + ending
            fl_kin = open(name_kin, 'r')
        TIMESTEPS = TIMESTEPS_HLRN
        hlrnflag = 1
    except:
        if d > NUM_OF_TRAJ_RZ:
            dummy=1
        else:
            name_kin = in_folder + str(jb) + '/' + str(d) + ending
            fl_kin = open(name_kin, 'r')
            TIMESTEPS = TIMESTEPS_RZ
            hlrnflag = 0

    return fl_kin, TIMESTEPS, hlrnflag


def Readfiles(ctr, d, jb, name, in_folder, TIMESTEPS_HLRN, TIMESTEPS_RZ, NUM_OF_TRAJ_RZ):
    flag = 0

    # open file for potential energy;
    # keep in mind that relevant information is only stored in every tenth line

    fl_pe, TIMESTEPS, hlrnflag = Openfile(d, in_folder, name, jb, ".lammpstrj", TIMESTEPS_RZ, TIMESTEPS_HLRN, NUM_OF_TRAJ_RZ)
    Traj = np.zeros([4, TIMESTEPS])
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

    fl_kin, TIMESTEPS, hlrnflag = Openfile(d, in_folder, name, jb, ".dat", TIMESTEPS_RZ, TIMESTEPS_HLRN, NUM_OF_TRAJ_RZ)

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
    fl_kin.close()
    fl_pe.close()
    return Traj, TIMESTEPS, hlrnflag
