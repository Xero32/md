####################################################################################################################
# Evaluate MD data
#
#   Track particle trajectories, calculate transition rates, compute analytical
#   solution to the Rate Equation Model, compute Residence Time
#
#   Calculate Energy Distribution, Angular Distribution
#
#   #NB#
'''How To:'''
#    python3 mdeval.py <Angle> <Temperature> <[--Options]>
#                   Options:    --hlrn
#                               --nrg (float)
#                               --start (int) [set start time for averaging of Transition Rates]
#                               --end (int) [set end time for averaging of Transition rates]
#
#   Program reads files from directories    /home/<user>/lammps/111/a{angle}t{temperature}e{energy}/{job}
#   or alternatively:                       /home/<user>/lammps/111/HLRN/a{angle}t{temperature}e{energy}/{job}
#   g.HOME = /home/<user> is defined in cfg.py, as are all important global variables
#
#   The program reads '.dat' and '.lammpstrj' files, where the file names themselves are the number (integer)
#   of trajectory in the corresponding job directory, in numerical order, beginning at '1.dat'
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
#       'mdeval':   contains general routines

#main_mdtraj.py
# -*- coding: utf-8 -*-
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import scipy, math, sys, time, argparse
# import plot, bounce, nrg, smooth, stats, params
# import cfg as g
# import matplotlib
import mdeval as m
#matplotlib.rcParams.update({'font.size': 16})
def printStd(a,b):
    print("%.4f +/- %.4f" %(a,b))

#
def main():
    m.ReadfileFn()
    m.StickingFn()
    m.PlotPopulations()
    TimePrime = m.Binning()
    m.PlotTransitionPopulations(TimePrime)
    TimePrime2 = m.Calc_BinTransitionRates()
    m.PlotTransitionRates(TimePrime2)
    m.Write_Const_TransitionRates()
    N_Matrix, Time2, multiplier = m.AnalyticalSolution()
    m.PlotSolution(N_Matrix, Time2, TimePrime, multiplier)
    m.Check_EquilibrationCondition(False)
    En_name = m.EnergyHandling()
    m.PlotEnergyDistr(En_name)
    bc = m.CalcTrajectoryDecay()
    m.Std_EnergyOverTime(bc)
    m.Plot_MeanEnergyOverTime()
    m.Find_EquilibrationTime()
    m.PlotTrajectoryDecay(bc)
    m.CalcAvgBounceTime()

if __name__ == "__main__":
    main()
