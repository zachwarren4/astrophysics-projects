import numpy as np
import pandas as pd
import math as math
import matplotlib.cbook as cbook
import matplotlib.pyplot as plt
import TreeCorr as tc

f, plts = plt.subplots(3,1,figsize=(15,15))


def correlation(data, rand,dataPoints):
    LocationD = r'C:\Users\Zachary Warren\OneDrive\S2 2015-2016\ASTR 3800\Homework3\src\%s' %data
    LocationR = r'C:\Users\Zachary Warren\OneDrive\S2 2015-2016\ASTR 3800\Homework3\src\%s' %rand

    dfD = pd.read_csv(LocationD, delimiter = ' ', names = ['RA', 'DEC', 'z'], skipinitialspace = True)
    dfR = pd.read_csv(LocationR, delimiter = ' ', names = ['RA', 'DEC', 'z'], skipinitialspace = True)



    bins = np.linspace(.1,20,15)

    corr=

    return (corr,bins)

