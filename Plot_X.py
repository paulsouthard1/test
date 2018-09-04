#%%
# coding: utf-8



import numpy as np
import plotly as py
import plotly.graph_objs as go
import os
import pickle
import scipy.special
import matplotlib.pyplot as plt
from scipy.stats import gamma

direc = os.path.dirname(os.path.realpath(__file__)) + "\\"



#Read in reach file to determine # of reaches
ReachFile = direc + "reach.txt"
with open(ReachFile,'r') as f:
    #Variable reachlist contains coordinates for each reach to be analyzed
    reachlist = f.readlines()

#Read in pickle file with number of sections
sect = pickle.load(open("sectnum","rb"))

#Blank dictionary
Elevs = {}
Dists = {}

# Load in profile data section-by-section and clip to valley walls
for i in range(len(reachlist)):
    sectnum = "mn_{}".format(i)
    for j in range(sect[sectnum]+1):
        # Read in profile data
        lookup = "CS_{}_{}".format(i,j+1) + ".txt"
        out_table = "T_" + lookup
        CS_Array = np.loadtxt(out_table,delimiter = ',',skiprows = 1,usecols = (1,2))
        # Store longitudinal distance and elevations as own variables
        Dist = CS_Array[:,0]
        Elev = CS_Array[:,1]
        # Working section: keep full profiles to view in plots
        Dist2 = Dist
        Elev2 = Elev
        # Go through each side of section and look for end of active corridor by several methods
        # Create arrays to loop through each side of section
        side = range(int((len(Elev)/2)+1))
        rside = np.flip(side,0)
        lside = np.arange(int(len(Elev)/2),(int(len(Elev))),1)
#        rslope = np.zeros(len(rside)+1)
#        lslope = np.zeros(len(lside)+1)
        rend = 0
        lend = len(Elev)
        #Find end of active corridor on right side
#        for k in rside[1:-1:]:
#            rslope[k] = np.absolute((Elev[k]-Elev[k-1])/(Dist[k]-Dist[k-1]))
        for k in rside[1:]:
            if Elev[k] >= (Elev[rside[0]] + 5):
                rend = k
                print("Total Elev")
                break
#            elif (rslope[k] >= 1.5) & (rslope[k-1] >= 1.5):
#                rend = k
#                print("Slope")
#                break
        #Find end of active corridor on left side
#        for k in lside[1:-1:]:
#            lslope[k-51] = np.absolute((Elev[k]-Elev[k-1])/(Dist[k]-Dist[k-1]))
        for k in lside[1:]:
            if Elev[k] >= (Elev[lside[0]] + 5):
                lend = k
                print("Total Elev")
                break
#            elif (lslope[(k-51)] >= 1.5) & (lslope[(k-51)+1] >= 1.5):
#                lend = k
#                print("Slope")
#                break
        savename = "set_{}".format(j)
        savename2 = "set2_{}".format(j)
        Elevs[savename] = Elev[rend:lend]
        Dists[savename] = Dist[rend:lend]
        Elevs[savename2] = Elev2
        Dists[savename2] = Dist2
for i in range(len(reachlist)):
    sectnum = "mn_{}".format(i)
    for j in range(sect[sectnum]+1):
        key = "set_{}".format(j)
        key2 = "set2_{}".format(j)
        trace1 = go.Scatter(
            x = Dists[key],
            y = Elevs[key],
            line = dict(
                color = ('rgb(255,69,0)'),
                width = 2)
        )
        trace2 = go.Scatter(
            x = Dists[key2],
            y = Elevs[key2],
            line = dict(
                color = ('rgb(0,128,128)'),
                width = 4)
        )
        data = [trace2,trace1]
        layout = go.Layout(
            autosize=False,
            width=2000,
            height=400,
            margin=go.Margin(
                l=100,
                r=50,
                b=100,
                t=100,
                pad=4
        ),
        title = 'Woodruff Canyon Channel Cross Sections {}'.format(j+1),
        font=dict(family='Arial', size=22),
        xaxis = dict(
            title = 'Distance Across Channel (m)',showline=True),
        yaxis = dict(
            title = 'Elevation Above Thalweg (m)',
            showline=True,
            range = [np.amin(Elevs[key2]),(np.amin(Elevs[key2])+10)]
            ),
        legend = dict(
            font = dict(
                size=19))
        )
        fig = go.Figure(data=data, layout=layout)
        filename = 'test_{}.html'.format(j)
        #py.offline.plot(fig,filename = filename)
#%%
#Save Data as pickle files
EFile = open(direc + "\\distribs\\efile","wb")
DFile = open(direc + "\\distribs\\dfile","wb")
pickle.dump(Elevs,EFile)
pickle.dump(Dists,DFile)

#%%
# Attempt to plot distributions and fit parameters


#Blank Dictionaries for Raw Height above thalwegs and normalized Height above thalweg
HATs = {}
S_HATs = {}

# Normaliz profile data X-sect by X-sect
for i in range(len(reachlist)):
    sectnum = "mn_{}".format(i)
    for j in range(sect[sectnum]+1):
        savename = "set_{}".format(j)
        Elev = Elevs[savename]
        Dist = Dists[savename]
        n = len(Elev)
        sampkey = "samp_{}_{}".format(i,j+1)
        Elevs[sampkey] = n
        thal = np.amin(Elev)
        HAT = Elev - thal
        S_HAT = HAT/(np.mean(HAT))
        HATs[savename] = HAT
        S_HATs[savename] = S_HAT
        
        
# Use Moving Window to get larger samples


# DEFINE WINDOW SIZES TO BE TESTED HERE, do odd #'s
Windows = [15]

WDists = {}

for i in range(len(Windows)):
    winsize = Windows[i]
    start = int(np.floor(winsize/2))
    end = sect[sectnum]-start
    extent = np.arange(start,end,1)
    for j in extent:
        winname = "{}_set_{}".format(winsize,j)
        S_winname = "S_{}_set_{}".format(winsize,j)
        window = np.arange((j-int(np.floor(winsize/2))),(j+int(np.floor(winsize/2)+1)),1)
        Ewin = []
        for k in window:
            savename = "set_{}".format(k)
            Ewin = np.concatenate((Ewin,HATs[savename]),axis = 0)
            S_Ewin = np.concatenate((Ewin,S_HATs[savename]),axis = 0)
        WDists[winname] = Ewin
        WDists[S_winname] = S_Ewin


#%% Some simple gamma dist tests
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.stats import gamma

#Creat function to plot freqdistribution values
def FreqDist(data,bins):
    counts = np.zeros(len(x))
    for i in range(len(x)):
        if i == 0:
            lower = 0
        else:
            lower= x[i-1]
        upper = x[i]
        for j in range(len(testdata)):
            if (testdata[j] > lower) & (testdata[j] <= upper):
                counts[i] = counts[i]+1
    freq = counts/np.sum(counts)
    return freq,counts

fig, ax = plt.subplots(1, 1)

testdata = WDists["15_set_10"]

fit_alpha, fit_loc, fit_beta=gamma.fit(testdata)
print(fit_alpha, fit_loc, fit_beta)

x = np.logspace(np.log10(0.01),np.log10(30))
data = gamma.pdf(x,2.8)
ax.plot(x, data, 'k-', lw=2, label='frozen pdf')

#Create frequency distribution with numpy
freq, counts = FreqDist(testdata,x)


ax.plot(x,freq,'r.')
ax.legend(loc='best', frameon=False)
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.001,1)
plt.show()

#%%
fig, ax = plt.subplots(1, 1)
a = 1.99
mean, var, skew, kurt = gamma.stats(a, moments='mvsk')


x = np.linspace(gamma.ppf(0.01, a),
                gamma.ppf(0.99, a), 100)
ax.plot(x, gamma.pdf(x, a),
       'r-', lw=5, alpha=0.6, label='gamma pdf')

rv = gamma(a)
ax.plot(x, rv.pdf(x), 'k-', lw=3, label='best-fit distribution')

#%%
testdata = WDists["15_set_10"]
x = np.logspace(np.log10(0.01),np.log10(10))

def FreqDist(data,bins):
    counts = np.zeros(len(bins))
    for i in range(len(bins)):
        if i == 0:
            lower = 0
        else:
            lower= bins[i-1]
        upper = bins[i]
        for j in range(len(data)):
            if (data[j] > lower) & (data[j] <= upper):
                counts[i] = counts[i]+1
    freq = counts/np.sum(counts)
    return freq,counts


FreqDist(testdata,x)

#%%

x = np.logspace(np.log10(0.01),np.log10(30))
alpha = np.linspace(0.4,10,97)

def MinError(freq,bins,alpha):
    from numpy import zeros,mean,argmin,absolute
    from scipy.stats import gamma
    MSE = zeros(len(alpha))
    for i in range(len(alpha)):
        a = alpha[i]
        SqErr = zeros(len(bins))
        GamDist = gamma.pdf(bins,a)
        for j in range(len(bins)):
            SqErr[j] = ((GamDist[j]-freq[j])/GamDist[j])**2
        print(SqErr)
        MSE[i] = mean(SqErr)
    ai = argmin(MSE)
    a_par = alpha[ai]
    return a_par,ai,MSE

a_par,ai,MSE = MinError(freq,x,alpha)


#%% Produce histograms with fitted gamma distributions for all windows

import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.stats import gamma

#%%
#Save Data as pickle files
EFile = open(direc + "\\distribs\\efile","wb")
DFile = open(direc + "\\distribs\\dfile","wb")
pickle.dump(Elevs,EFile)
pickle.dump(Dists,DFile)

#Create function to plot frequency distribution of values
def FreqDist(data,bins):
    counts = np.zeros(len(x))
    for i in range(len(x)):
        if i == 0:
            lower = 0
        else:
            lower= x[i-1]
        upper = x[i]
        for j in range(len(testdata)):
            if (testdata[j] > lower) & (testdata[j] <= upper):
                counts[i] = counts[i]+1
    freq = counts/np.sum(counts)
    return freq,counts

#Create function to find best fit for frequency distribution
def MinError(freq,bins,alpha):
    from numpy import zeros,mean,argmin,absolute
    from scipy.stats import gamma
    MSE = zeros(len(alpha))
    for i in range(len(alpha)):
        a = alpha[i]
        SqErr = zeros(len(bins))
        GamDist = gamma.pdf(bins,a)
        for j in range(len(bins)):
            SqErr[j] = ((GamDist[j]-freq[j])/GamDist[j])**2
        print(SqErr)
        MSE[i] = mean(SqErr)
    ai = argmin(MSE)
    a_par = alpha[ai]
    return a_par,ai,MSE

#Create loop to produce plots for all 
Windows = [9,15,25]

WDists = {}
sectnum = 26

for i in range(len(Windows)):
    winsize = Windows[i]
    start = int(np.floor(winsize/2))
    end = sectnum-start
    extent = np.arange(start,end,1)
    for j in extent:
        winname = "{}_set_{}".format(winsize,j)
        S_winname = "S_{}_set_{}".format(winsize,j)
        window = np.arange((j-int(np.floor(winsize/2))),(j+int(np.floor(winsize/2)+1)),1)
        Ewin = []
        for k in window:
            savename = "set_{}".format(k)
            Ewin = np.concatenate((Ewin,HATs[savename]),axis = 0)
            S_Ewin = np.concatenate((Ewin,S_HATs[savename]),axis = 0)
        WDists[winname] = Ewin
        WDists[S_winname] = S_Ewin




fig, ax = plt.subplots(1, 1)
ax.plot(x,freq,'r.')
ax.legend(loc='best', frameon=False)
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.001,1)
plt.show()
