import numpy as np
import pandas as pd
import math as math
import matplotlib.cbook as cbook
import matplotlib.pyplot as plt
from scipy import spatial

f, plts = plt.subplots(3,1,figsize=(15,15))

def toX(row):
    dist = 3000*row['z']
    return dist*math.cos(math.radians(row['DEC']))*math.cos(math.radians(row['RA']))

def toY(row):
    dist = 3000*row['z']
    return dist*math.cos(math.radians(row['DEC']))*math.sin(math.radians(row['RA']))

def toZ(row):
    dist = 3000*row['z']
    return dist*math.sin(math.radians(row['DEC']))

def dist(loc1,loc2,df):
    return np.sqrt((df['xVal'].iloc[loc1] - df['xVal'].iloc[loc2])**2 + (df['yVal'].iloc[loc1] - df['yVal'].iloc[loc2])**2 + (df['zVal'].iloc[loc1] - df['zVal'].iloc[loc2])**2)

def countForBin(bin,df,dataPoints):
    finalCount = np.linspace(0,0,15)
    y = 0
    while y < dataPoints + 2:
        for x in range(y + 1,dataPoints):
            w = int(round((dist(y,x,df)+1)/1.436))
            if w < 15:
                finalCount[w]+=1
        y = y + 1
    return finalCount

def getCorrelation(dCount,rCount,nD,rD):
        ratio = [float(dCount)/rCount for dCount,rCount in zip(dCount,rCount)]
        ratio[:] = [(((nD/rD)**2)* x) - 1 for x in ratio]
        return ratio

def correlation(data, rand,dataPoints):
    LocationD = r'C:\Users\Zachary Warren\OneDrive\S2 2015-2016\ASTR 3800\Homework3\src\%s' %data
    LocationR = r'C:\Users\Zachary Warren\OneDrive\S2 2015-2016\ASTR 3800\Homework3\src\%s' %rand

    dfD = pd.read_csv(LocationD, delimiter = ' ', names = ['RA', 'DEC', 'z'], skipinitialspace = True)
    dfR = pd.read_csv(LocationR, delimiter = ' ', names = ['RA', 'DEC', 'z'], skipinitialspace = True)

    dfD['xVal'] = dfD.apply(lambda row: toX(row), axis = 1)
    dfD['yVal'] = dfD.apply(lambda row: toY(row), axis = 1)
    dfD['zVal'] = dfD.apply(lambda row: toZ(row), axis = 1)

    dfR['xVal'] = dfR.apply(lambda row: toX(row), axis = 1)
    dfR['yVal'] = dfR.apply(lambda row: toY(row), axis = 1)
    dfR['zVal'] = dfR.apply(lambda row: toZ(row), axis = 1)

    dfD.to_csv('data4.csv')

    dfDs = dfD.sample(5000)
    dfRs = dfR.sample(5000)
    dfDs = dfDs.reset_index()
    dfRs = dfRs.reset_index()

    bins = np.linspace(.1,20,15)
    print(bins)
    countD = countForBin(bins,dfDs,dataPoints)
    countR = countForBin(bins,dfRs,dataPoints)

    print('Count is ',countD, 'countR ', countR)

    corr= getCorrelation(countD,countR,dataPoints,dataPoints)
    print(corr)

    return (corr,bins)


def correlationDM(data, rand,dataPoints):
    LocationD = r'C:\Users\Zachary Warren\OneDrive\S2 2015-2016\ASTR 3800\Homework3\src\%s' %data
    LocationR = r'C:\Users\Zachary Warren\OneDrive\S2 2015-2016\ASTR 3800\Homework3\src\%s' %rand

    dfD = pd.read_csv(LocationD, delimiter = ' ', names = ['xVal', 'yVal', 'zVal'], skipinitialspace = True)
    dfR = pd.read_csv(LocationR, delimiter = ' ', names = ['xVal', 'yVal', 'zVal'], skipinitialspace = True)

    dfDs = dfD.sample(5000)
    dfRs = dfR.sample(5000)

    dfDs = dfDs.reset_index()
    dfRs = dfRs.reset_index()

    bins = np.linspace(.1,20,15)
    print(bins)
    countD = countForBin(bins,dfDs,dataPoints)
    countR = countForBin(bins,dfRs,dataPoints)

    print('Count is ',countD, 'countR ', countR)

    corr= getCorrelation(countD,countR,dataPoints,dataPoints)
    print(corr)

    return (corr,bins)


def getBias(dCount,rCount):
        rCount[:] = [x + 1 for x in rCount]
        ratio = [float(dCount)/rCount for dCount,rCount in zip(dCount,rCount)]
        ratio[:] = [np.sqrt(x) for x in ratio]
        return ratio

corr20,bins20 = correlation('SDSS_Mr20_rspace.dat','SDSS_random.dat',5000)
corr21,bins21 = correlation('SDSS_Mr21_rspace.dat','SDSS_random.dat',5000)
corr20z,bins20z = correlation('SDSS_Mr20_zspace.dat','SDSS_random.dat',5000)
corrDM,binsDM = correlationDM('DM.dat','DM_random.dat',5000)
bias20 = getBias(corr20,corrDM)
bias21 = getBias(corr21,corrDM)


plts[0].loglog(bins20, corr20)
plts[0].loglog(bins21,corr21, 'r')
plts[1].loglog(bins20,corr20)
plts[1].loglog(bins20z,corr20z,'r')
plts[2].plot(np.log10(bins20),bias20)
plts[2].plot(np.log10(bins21),bias21,'r')

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