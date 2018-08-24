import numpy as np
import math, sys
import matplotlib.pyplot as plt
import plot
import smooth
import compareAng as ca

Theta = np.full((6, 10000), -9999.)

######################################################################

ang = 30
temp= 300
bnc = 4

ca.Assertion(ang, temp, bnc)

ca.TempCompare(ang,bnc, Theta)
ca.AngleCompare(temp,bnc, Theta)
ca.BounceCompare(temp,ang, Theta)
