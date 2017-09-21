import numpy as np
import math as math
import matplotlib.pyplot as plt
from scipy import spatial

f, plts = plt.subplots(3,1,figsize=(15,15))

#constant bins for easy adjustment
BINSnotlog = np.linspace(-1,1.3,16)
BINS = np.power(10,BINSnotlog)

#create X,Y,Z coordinates
def toXYZ(RA,DEC,z):
    x = []
    y = []
    zout = []
    for i in range(0,len(RA)):
        r = 3000*z[i]
        x.append(r*math.cos(math.radians(DEC[i]))*math.cos(math.radians(RA[i])))
        y.append(r*math.cos(math.radians(DEC[i]))*math.sin(math.radians(RA[i])))
        zout.append(r*math.sin(math.radians(DEC[i])))
    return list(zip(x,y,zout))

#get the bias
def getBias(corr,corrDM):
    bias =[np.sqrt(float(corr)/corrDM) for corr,corrDM in zip(corr,corrDM)]
    return bias

#calculate the correlation, change to x,y,z coordinates only if using non-DM data
#uses kdTrees to count neighbors
def correlation(data, rand, DM, bins):
    if (DM==False):
        RA,DEC,z = np.loadtxt(data, unpack=True)
        RAR,DECR,zr = np.loadtxt(rand, unpack=True)

        xyz = toXYZ(RA,DEC,z)
        xyzr = toXYZ(RAR,DECR,zr)

        #length of list
        numD = len(RA)
        numR = len(RAR)

    else:
        x,y,z = np.loadtxt(data,unpack=True)
        xr,yr,zr = np.loadtxt(rand,unpack=True)

        xyz=list(zip(x,y,z))
        xyzr=list(zip(xr,yr,zr))

        #length of list
        numD = len(x)
        numR = len(xr)

    kdD = spatial.cKDTree(xyz)
    kdR = spatial.cKDTree(xyzr)

    countD=[]
    countR=[]

    #normalizing the correlation function
    norm = ((numR/numD)**2)

    #append counts for each bin from r=0 to r=x
    for x in bins:
        countD.append(kdD.count_neighbors(kdD,x))
        countR.append(kdR.count_neighbors(kdR,x))

    # subtracting out length of list to eliminate the x1,x1
    #comparison that says all points are neighbors with themselves
    countD[:] = [x-numD for x in countD]
    countR[:] = [x-numR for x in countR]

    #subtracts all counts from bins less than current bin
    for i in range(1,16):
        countD[i] = countD[i] - sum(countD[0:i])
        countR[i] = countR[i] - sum(countR[0:i])

    #calculates correlation function
    corr = [(norm*(float(countD)/(countR+1))-1) for countD,countR in zip(countD,countR)]

    return (corr)

#get correlations and bias
corr20 = correlation('SDSS_Mr20_rspace.dat', 'SDSS_random.dat',DM=False, bins=BINS)
corr21 = correlation('SDSS_Mr21_rspace.dat', 'SDSS_random.dat',DM=False, bins=BINS)
corr20z = correlation('SDSS_Mr20_zspace.dat', 'SDSS_random.dat',DM=False, bins=BINS)
corrDM = correlation('DM.dat','DM_random.dat', DM=True, bins=BINS)
bias20 = getBias(corr20, corrDM)
bias21 = getBias(corr21, corrDM)

#set up plots
plts[0].plot(np.log10(BINS), np.log10(corr20))
plts[0].plot(np.log10(BINS), np.log10(corr21), 'r')
plts[1].plot(np.log10(BINS), np.log10(corr20))
plts[1].plot(np.log10(BINS), np.log10(corr20z),'r')
plts[2].plot(np.log10(BINS),bias20)
plts[2].plot(np.log10(BINS),bias21,'r')

plts[0].set_title('Correlation vs Distance')
plts[0].set_xlabel('log(r)')
plts[0].set_ylabel('log(E(r))')
plts[0].legend(['Mr20','Mr21'])

plts[1].set_title('Correlation vs Distance')
plts[1].set_xlabel('log(r)')
plts[1].set_ylabel('log(E(r))')
plts[1].legend(['Mr20 Real','Mr20 Z'])

plts[2].set_xlabel('log(r)')
plts[2].set_ylabel('b(r)')
plts[2].set_title('Bias vs. Distance')
plts[2].legend(['Mr20','Mr21'])

plt.show()


