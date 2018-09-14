# input.py
import numpy as np
import bounce
import argparse
import cfg as g

def InpParams(parser):
    parser.add_argument("a", help='Angle', default='30', nargs='?')
    parser.add_argument("t", help='Temperature', default='300', nargs='?')
    parser.add_argument("--fp", help='Temporal_Fixpoint for initial value problem')
    parser.add_argument("--enrg", help='Incident Energy')
    parser.add_argument("--start", help="Sets start time to calculate average transition rate", default=20)
    parser.add_argument("--end", help="Sets end time to calculate average transition rate", default=30)
    parser.add_argument("--hlrn", help="Use trajectories from HLRN runs", action='store_true')
    args = parser.parse_args()
    if args.fp:
        fp = float(args.fp)
        print(fp)
    else:
        fp = int(20//0.025) #default for t=20ps at 0.00025 ps timestep with data written out every 100 steps
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
    if args.enrg:
        energy = args.enrg
    else:
        energy = args.a
    if args.hlrn:
        g.HLRN = 1
    else:
        g.HLRN = 0
    return startps, endps, angle, temperature, energy, g.HLRN, fp

def Parameters(temperature, energy): #, HLRN):
    if int(temperature) == 80:
        if g.HLRN:
            g.TIMESTEPS_RZ = 120000 // 100
            g.TIMESTEPS_HLRN = 160000 // 100
            g.NUM_OF_TRAJ_RZ = 200
            g.NUM_OF_TRAJ = 1000
            g.nbnc = 40
            if float(energy) == 25.90:
                g.jobs = (1998313, 0)
                g.nbnc = 30
                g.TIMESTEPS_HLRN = 120000 // 100
            elif float(energy) == 51.80:
                g.jobs = (1998314, 0)
                g.nbnc = 30
                g.TIMESTEPS_HLRN = 120000 // 100
            elif float(energy) == 103.60:
                g.jobs = (1998315, 0)
                g.nbnc = 30
                g.TIMESTEPS_HLRN = 120000 // 100
            else:
                g.jobs = (1931936, 1931937)
        else:
            g.TIMESTEPS_RZ = 120000 // 100
            g.TIMESTEPS_HLRN = 120000 // 100
            g.NUM_OF_TRAJ_RZ = 200
            g.NUM_OF_TRAJ = 200
            g.nbnc = 30
            g.jobs = (7094790, 7107184, 7107185, 7107186, 7107187, 7107195, 7131096, 7131097, 7131115, 7131116)
        #jobs = (7279669, 7279673, 7279675, 7279676, 7279677)
        te = 11 // 0.025 #T=80 K
    elif int(temperature) == 120:
        if g.HLRN:
            g.TIMESTEPS_RZ = 120000 // 100
            g.TIMESTEPS_HLRN = 160000 // 100
            g.NUM_OF_TRAJ_RZ = 200
            g.NUM_OF_TRAJ = 1000
            g.nbnc = 40
            g.jobs = (1931939, 1935306)
        else:
            g.TIMESTEPS_RZ = 120000 // 100
            g.TIMESTEPS_HLRN = 120000 // 100
            g.NUM_OF_TRAJ_RZ = 200
            g.NUM_OF_TRAJ = 200
            g.nbnc = 30
            g.jobs = ()
    elif int(temperature) == 160:
        if g.HLRN:
            g.TIMESTEPS_RZ = 120000 // 100
            g.TIMESTEPS_HLRN = 160000 // 100
            g.NUM_OF_TRAJ_RZ = 200
            g.NUM_OF_TRAJ = 1000
            g.nbnc = 40
            g.jobs = (1931943, 1935310)

        else:
            g.TIMESTEPS_RZ = 120000 // 100
            g.TIMESTEPS_HLRN = 120000 // 100
            g.NUM_OF_TRAJ_RZ = 200
            g.NUM_OF_TRAJ = 200
            g.nbnc = 30
            g.jobs = ()
    elif int(temperature) == 190:
        if g.HLRN:
            g.TIMESTEPS_RZ = 120000 // 100
            g.TIMESTEPS_HLRN = 160000 // 100
            g.NUM_OF_TRAJ_RZ = 200
            g.NUM_OF_TRAJ = 1000
            g.nbnc = 40
            if float(energy) == 25.90:
                g.jobs = (1998316, 0)
                g.nbnc = 30
                g.TIMESTEPS_HLRN = 120000 // 100
            elif float(energy) == 51.80:
                g.jobs = (1998317, 0)
                g.nbnc = 30
                g.TIMESTEPS_HLRN = 120000 // 100
            elif float(energy) == 103.60:
                g.jobs = (1998318, 0)
                g.nbnc = 30
                g.TIMESTEPS_HLRN = 120000 // 100
            else:
                g.jobs = (1935315, 1935316)
        else:
            g.TIMESTEPS_RZ = 120000 // 100
            g.TIMESTEPS_HLRN = 120000 // 100
            g.NUM_OF_TRAJ_RZ = 200
            g.NUM_OF_TRAJ = 200
            g.nbnc = 30
            g.jobs = (7237382, 7237376, 7237379, 7237380, 7237381, 8199182, 8199181, 8199180, 8199179, 8199176)
    elif int(temperature) == 240:
        if g.HLRN:
            g.TIMESTEPS_RZ = 120000 // 100
            g.TIMESTEPS_HLRN = 160000 // 100
            g.NUM_OF_TRAJ_RZ = 200
            g.NUM_OF_TRAJ = 1000
            g.nbnc = 40
            g.jobs = (1931944, 1931945)
        else:
            g.TIMESTEPS_RZ = 200000 // 100
            g.TIMESTEPS_HLRN = 200000 // 100
            g.NUM_OF_TRAJ_RZ = 1000
            g.NUM_OF_TRAJ = 1000
            g.nbnc = 50
            g.jobs = (1967541, 1967542, 1969467, 1969466)
    elif int(temperature) == 300:
        g.TIMESTEPS_RZ = 120000 // 100
        g.NUM_OF_TRAJ_RZ = 200
        BOTH = 0
        if g.HLRN == 0:
            g.TIMESTEPS_HLRN = 120000 // 100
            g.NUM_OF_TRAJ = 200
            g.nbnc = 30
            if float(energy) != 70.0:
                #TODO remove clip in job array
                g.jobs = (7094785, 7107197, 7107198, 7107199, 7107200, 7107201, 7131101, 7131102, 7131113, 7131114)
            else:
                g.jobs = (7227490, 7227494, 7227495, 7227496, 7227497)
        else:
            g.TIMESTEPS_HLRN = 200000 // 100
            g.NUM_OF_TRAJ = 500
            g.nbnc = 50
            if BOTH == 1:
                g.jobs = (7094785, 7107197, 7107198, 7107199, 7107200, 7107201, 7131101, 7131102, 7131113, 7131114, 1922961, 1922962, 1922963, 1922964) # rz + hlrn jobs
            else:
                if float(energy) == 25.90:
                    g.TIMESTEPS_HLRN = 120000 // 100
                    g.NUM_OF_TRAJ = 1000
                    g.jobs = (1998319, 0)
                    g.nbnc = 30
                elif float(energy) == 51.80:
                    g.TIMESTEPS_HLRN = 120000 // 100
                    g.NUM_OF_TRAJ = 1000
                    g.jobs = (1998320, 0)
                    g.nbnc = 30
                elif float(energy) == 103.60:
                    g.TIMESTEPS_HLRN = 120000 // 100
                    g.NUM_OF_TRAJ = 1000
                    g.jobs = (1998321, 0)
                    g.nbnc = 30
                else:
                    g.jobs = (1922961, 1922962, 1922963, 1922964, 1965997, 1966017, 1965996, 1966006)

 # hlrn jobs

    #return TIMESTEPS_RZ, TIMESTEPS_HLRN, NUM_OF_TRAJ_RZ, NUM_OF_TRAJ, nbnc, jobs

def Openfile(d, in_folder, in_folder_hlrn, name, jb, ending): #, TIMESTEPS_RZ, TIMESTEPS_HLRN, NUM_OF_TRAJ_RZ):
    try:
        try:
            name_kin = in_folder + str(jb) + '.bbatch.hsn.hlrn.de/' + str(d) + ending
            fl_kin = open(name_kin, 'r')
        except:
            name_kin = in_folder_hlrn + name + '/' + str(jb) + '.bbatch.hsn.hlrn.de/' + str(d) + ending
            fl_kin = open(name_kin, 'r')
        g.TIMESTEPS = g.TIMESTEPS_HLRN
        hlrnflag = 1
    except:
        if d > g.NUM_OF_TRAJ_RZ:
            dummy=1
        else:
            name_kin = in_folder + str(jb) + '/' + str(d) + ending
            fl_kin = open(name_kin, 'r')
            g.TIMESTEPS = g.TIMESTEPS_RZ
            hlrnflag = 0

    return fl_kin, g.TIMESTEPS, hlrnflag


def Readfiles(ctr, d, jb, name, in_folder, in_folder_hlrn): #, TIMESTEPS_HLRN, TIMESTEPS_RZ, NUM_OF_TRAJ_RZ):
    flag = 0

    # open file for potential energy;
    # keep in mind that relevant information is only stored in every tenth line

    fl_pe, g.TIMESTEPS, hlrnflag = Openfile(d, in_folder, in_folder_hlrn, name, jb, ".lammpstrj") #, TIMESTEPS_RZ, TIMESTEPS_HLRN, NUM_OF_TRAJ_RZ)
    Traj = np.zeros([7, g.TIMESTEPS])
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

    fl_kin, g.TIMESTEPS, hlrnflag = Openfile(d, in_folder, in_folder_hlrn, name, jb, ".dat") #, TIMESTEPS_RZ, TIMESTEPS_HLRN, NUM_OF_TRAJ_RZ)

    for a,line in enumerate(fl_kin):
        b = a - 1
        if a == 0:
            continue

        t, x, y, Traj[6,b], Traj[4,b], Traj[5,b], Traj[3,b], en = line.split()
        Traj[0,b] = bounce.normalKE(float(Traj[3,b]))
        Traj[1,b] = bounce.parallelKE(float(Traj[4,b]), float(Traj[5,b]))
        assert(Traj[0,b] >= 0)
        assert(Traj[1,b] >= 0)
    #end for (kinetic file)
    fl_kin.close()
    fl_pe.close()
    return Traj, g.TIMESTEPS, hlrnflag
