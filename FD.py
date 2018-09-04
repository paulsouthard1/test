# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 13:59:52 2018

@author: ps29626
"""

import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
from scipy.stats import gamma
import seaborn as sns



dir_path = os.path.dirname(os.path.realpath(__file__)) + "\\"

#Find name of containing folder, i.e. sample name
beg = dir_path.rfind("\\",0,len(dir_path)-2) +1
end = len(dir_path)-1
subsamp = dir_path[beg:end]




#USER SET FOLDER NAME HERE - Name your analysis 
foldername = "Initial"
analysispath = dir_path + foldername + "\\"

if not os.path.exists(analysispath):
    os.makedirs(analysispath)

#%%
#Define Modules

#Create function to plot frequency distribution of values
def FreqDist(data,bins):
    from numpy import zeros,sum
    counts = zeros(len(bins))
    for i in range(len(bins)):
        if i == 0:
            lower = 0
        else:
            lower= bins[i-1]
        upper = bins[i]
        for j in range(len(data)):
            if (data[j] > lower) & (data[j] <= upper):
                counts[i] = counts[i]+1
    freq = counts/sum(counts)
    return freq,counts

#Create function to select active corridor portion of each cross-section
def ActiveCorr(Dist,Elev,Height):
    from numpy import flip, arange
    Dist_U = Dist
    Elev_U = Elev
    side = range(int((len(Elev)/2)+1))
    rside = flip(side,0)
    lside = arange(int(len(Elev)/2),(int(len(Elev))),1)
    rend = 0
    lend = len(Elev)
    for k in rside[1:]:
        if Elev[k] >= (Elev[rside[0]] + Height):
            rend = k
            break
    for k in lside[1:]:
        if Elev[k] >= (Elev[lside[0]] + Height):
            lend = k
            break
    Dist_P = Dist[rend:lend]
    Elev_P = Elev[rend:lend]
    return Dist_U,Elev_U,Dist_P,Elev_P

#Create function to calculate dimensionless bed elevations for each cross-section
def DimBed(Elev):
    from numpy import mean, amin
    n = len(Elev)
    thal = amin(Elev)
    HAT = Elev - thal
    avg = mean(HAT)
    S_HAT = HAT/avg
    return HAT, S_HAT, n, avg

#Create function to use moving windows to analyze groups of multiple cross-sections
def Win(Sizes,ElevDict,WinDict,NumSections):
    from numpy import ceil, floor, arange
    for w in range(len(Sizes)):
        winsize = Sizes[w]
        number = np.ceil(NumSections/winsize)
        start = int(np.floor(winsize/2))
        end = (NumSections-start)
        extent = np.arange(start,end,1)
        for e in extent:
            winname = "{}_{}_win_{}".format(i,winsize,e)
            S_winname = "S_{}_{}_win_{}".format(i,winsize,e)
            window = np.arange((e-int(np.floor(winsize/2))),(e+int(np.floor(winsize/2)+1)),1)
            Ewin = []
            S_Ewin = []
            for z in window:
                lookup = "set_{}_{}".format(i,z)
                lookup_S = lookup + "_S"
                Ewin = np.concatenate((Ewin,ElevDict[lookup]),axis = 0)
                S_Ewin = np.concatenate((S_Ewin,ElevDict[lookup_S]),axis = 0)
            WinDict[winname] = Ewin
            WinDict[S_winname] = S_Ewin
    return WinDict

#Create function to find best fit for frequency distribution
def MinError(freq,bins,alpha,loc,scale):
    from numpy import zeros,mean,argmin,absolute
    from scipy.stats import gamma
    MSE = zeros(len(alpha))
    for i in range(len(alpha)):
        a = alpha[i]
        SqErr = zeros(len(bins))
        GamDist = gamma.pdf(bins,a,loc = loc,scale = scale)
        for j in range(len(bins)):
            if freq[j] == 0:
                SqErr[j] == 0
            else:
                SqErr[j] = ((GamDist[j]-freq[j])/freq[j])**2
        MSE[i] = np.sum(SqErr)
    ai = argmin(MSE)
    a_par = alpha[ai]
    return a_par,ai,MSE

#%%
#Read in reach file to determine # of reaches
ReachFile = dir_path + "reach_" + subsamp + ".txt"
with open(ReachFile,'r') as f:
    #Variable reachlist contains coordinates for each reach to be analyzed
    reachlist = f.readlines()

#Read in pickle file with number of sections
sect = pickle.load(open("sectnum","rb"))

#Create Blank Dictionaries to store raw cross-section datasets and active corridor datasets
Elevs = {}
Dists = {}

# Load in profile data section-by-section and clip to valley walls
for i in range(len(reachlist)):
    sectnum = "mn_{}".format(i)
    for j in range(sect[sectnum]+1):
        #USER SET ACTIVE CORRIDOR HEIGHT HERE
        Height = 5
        # Read in profile data
        lookup = "CS_{}_{}".format(i,j+1) + ".txt"
        out_table = "T_" + lookup
        CS_Array = np.loadtxt(out_table,delimiter = ',',skiprows = 1,usecols = (1,2))
        # Store longitudinal distance and elevations as own variables
        Dist = CS_Array[:,0]
        Elev = CS_Array[:,1]
        Dist_U,Elev_U,Dist_P,Elev_P = ActiveCorr(Dist,Elev,Height)
        savename = "set_{}_{}".format(i,j)
        savename_U = "set_U_{}_{}".format(i,j)
        Dists[savename] = Dist_P
        Dists[savename_U] = Dist_U
        Elevs[savename] = Elev_P
        Elevs[savename_U] = Elev_U

#%%
#Change Directory for storage of analysis files
os.chdir(analysispath)

#Dump Raw Cross-section and Active Corridor datasets out to pickle files
file = open(analysispath + "Efile", "wb")
pickle.dump(Elevs,file)
file.close()

file = open(analysispath + "Dfile", "wb")
pickle.dump(Dists,file)
file.close()

#Create Folder to store plots of Active Corridor Selection
foldername = "ACorr"
plotpath = analysispath + foldername + "\\"

if not os.path.exists(plotpath):
    os.makedirs(plotpath)
    
#Change Directory for storage of plots
os.chdir(plotpath)


#USER SET PLOTTING PREFERENCES HERE
Plot = False
#Plot selected active corridor on top of raw cross-section elevations
for i in range(len(reachlist)):
    if Plot == True:
        sectnum = "mn_{}".format(i)
        for j in range(sect[sectnum]+1):
            savename = "set_{}_{}".format(i,j)
            savename_U = "set_U_{}_{}".format(i,j)
            plt.style.use('seaborn-white')
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['font.serif'] = 'Ubuntu'
            plt.rcParams['font.monospace'] = 'Ubuntu Mono'
            plt.rcParams['font.size'] = 10
            plt.rcParams['axes.labelsize'] = 10
            plt.rcParams['axes.labelweight'] = 'bold'
            plt.rcParams['axes.titlesize'] = 10
            plt.rcParams['xtick.labelsize'] = 8
            plt.rcParams['ytick.labelsize'] = 8
            plt.rcParams['legend.fontsize'] = 10
            plt.rcParams['figure.titlesize'] = 12
            width, height = plt.figaspect(0.5)
            fig = plt.figure(figsize=(width,height), dpi=400)
            fig, ax = plt.subplots(1, 1)
            ax.plot(Dists[savename_U], Elevs[savename_U], linewidth=4, linestyle='-', label='Cross-section')
            ax.plot(Dists[savename], Elevs[savename], linewidth=3, linestyle='--', label='Active Corridor')
            ax.legend(frameon=False)
            plt.xlabel("Distance")
            plt.ylabel("Elevation")
            plt.ylim(np.amin(Elevs[savename]),np.amin(Elevs[savename])+15)
            filename = 'AC_{}_{}'.format(i,j+1)
            plt.title(filename)
            plt.savefig(filename)
            plt.close(fig)
    else:
        break

#%%
#Create blank dictionary for Dimensionless Elevations
HATs = {}

# Calculate height above thalweg and dimensionless bed elevation
for i in range(len(reachlist)):
    sectnum = "mn_{}".format(i)
    for j in range(sect[sectnum]+1):
        savename = "set_{}_{}".format(i,j)
        savename_S = savename + "_S"
        nsave = "n_{}_{}".format(i,j)
        avgsave = "avg_{}_{}".format(i,j)
        Elev = Elevs[savename]
        HAT, S_HAT, n, avg = DimBed(Elev)
        HATs[savename] = HAT
        HATs[savename_S] = S_HAT
        HATs[nsave] = n
        HATs[avgsave] = avg

#Change Directory for storage of analysis files
os.chdir(analysispath)

#Dump Dimensionless Elevation datasets out to pickle file
file = open(analysispath + "HATfile", "wb")
pickle.dump(HATs,file)
file.close()


#%%

#Create Frequency Distributions for individual cross-sections

#Create blank dictionary to store individual FD's
iFD = {}
iFDc = {}

#Create bins for Freq Dist storage
bins = np.linspace(0,5,num = 51)
bins_S = np.logspace(np.log10(0.01),np.log10(30))

for i in range(len(reachlist)):
    sectnum = "mn_{}".format(i)
    for j in range(sect[sectnum]+1):
        savename = "set_{}_{}".format(i,j)
        savename_S = savename + "_S"
        HAT = HATs[savename]
        SHAT = HATs[savename_S]
        freq,counts = FreqDist(HAT,bins)
        freq_S,counts_S = FreqDist(SHAT,bins_S)
        iFD[savename] = freq
        iFD[savename_S] = freq_S
        iFDc[savename] = counts
        iFDc[savename_S] = counts_S

#Dump Frequency and Count datasets out to pickle files
file = open(analysispath + "iFDFile", "wb")
pickle.dump(iFD,file)
file.close()

file = open(analysispath + "iFDcFile", "wb")
pickle.dump(iFDc,file)
file.close()

#Create Folder to store plots of 
foldername = "iFD"
plotpath = analysispath + foldername + "\\"

if not os.path.exists(plotpath):
    os.makedirs(plotpath)
    
#Change Directory for storage of plots
os.chdir(plotpath)


#USER SET PLOTTING PREFERENCES HERE
Plot = False
#Plot frequency distributions with HAT
for i in range(len(reachlist)):
    if Plot == True:
        sectnum = "mn_{}".format(i)
        for j in range(sect[sectnum]+1):
            savename = "set_{}_{}".format(i,j)
            plt.style.use('seaborn-white')
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['font.serif'] = 'Ubuntu'
            plt.rcParams['font.monospace'] = 'Ubuntu Mono'
            plt.rcParams['font.size'] = 10
            plt.rcParams['axes.labelsize'] = 10
            plt.rcParams['axes.labelweight'] = 'bold'
            plt.rcParams['axes.titlesize'] = 10
            plt.rcParams['xtick.labelsize'] = 8
            plt.rcParams['ytick.labelsize'] = 8
            plt.rcParams['legend.fontsize'] = 10
            plt.rcParams['figure.titlesize'] = 12
            width, height = plt.figaspect(0.5)
            fig = plt.figure(figsize=(width,height), dpi=400)
            fig, ax = plt.subplots(1, 1)
            ax.plot(bins, iFD[savename], 'bo', label='Cross-section')
            plt.xlabel("HAT")
            plt.ylabel("Frequency")
            plt.ylim(0,0.2)
            filename = 'FD_{}_{}'.format(i,j+1)
            plt.title(filename + " for HAT")
            plt.savefig(filename)
            plt.close(fig)
    else:
        break
    
#USER SET PLOTTING PREFERENCES HERE
Plot = False
#Plot frequency distributions with S_HAT's
for i in range(len(reachlist)):
    if Plot == True:
        sectnum = "mn_{}".format(i)
        for j in range(sect[sectnum]+1):
            savename = "set_{}_{}_S".format(i,j)
            plt.style.use('seaborn-white')
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['font.serif'] = 'Ubuntu'
            plt.rcParams['font.monospace'] = 'Ubuntu Mono'
            plt.rcParams['font.size'] = 10
            plt.rcParams['axes.labelsize'] = 10
            plt.rcParams['axes.labelweight'] = 'bold'
            plt.rcParams['axes.titlesize'] = 10
            plt.rcParams['xtick.labelsize'] = 8
            plt.rcParams['ytick.labelsize'] = 8
            plt.rcParams['legend.fontsize'] = 10
            plt.rcParams['figure.titlesize'] = 12
            width, height = plt.figaspect(0.5)
            fig = plt.figure(figsize=(width,height), dpi=400)
            fig, ax = plt.subplots(1, 1)
            ydata = iFD[savename]
            ax.plot(bins_S, ydata, 'bo', label='Cross-section')
            plt.xlabel("S_HAT")
            plt.ylabel("Frequency")
            plt.ylim(np.amin(ydata),np.amax(ydata))
            plt.xlim(np.amin(bins_S),np.amax(bins_S))
            plt.xscale('log')
            plt.yscale('log')
            filename = 'FD_{}_{}_S'.format(i,j+1)
            plt.title(filename + " for SHAT")
            plt.savefig(filename)
            plt.close(fig)
    else:
        break
    
#%%  Analysis of data via moving windows

#USER SET WINDOW SIZE HERE
Sizes = [5,15,25]
#Create Blank dictionary to store window samples
WD = {}
  
for i in range(len(reachlist)):
    sectnum = "mn_{}".format(i)
    NumSections = sect[sectnum]
    WinDict = Win(Sizes,HATs,WD,NumSections)

file = open(analysispath + "WFile", "wb")
pickle.dump(WD,file)
file.close()

#Create blank dictionary to store window FD's

wFD = {}
wFDc = {}

#Run loop to perform frequency distribution for window calculations
for i in range(len(reachlist)):
    sectnum = "mn_{}".format(i)
    NumSections = sect[sectnum]
    for w in range(len(Sizes)):
        winsize = Sizes[w]
        number = np.ceil(NumSections/winsize)
        start = int(np.floor(winsize/2))
        end = (NumSections-start)
        extent = np.arange(start,end,1)
        for e in extent:
            winname = "{}_{}_win_{}".format(i,winsize,e)
            S_winname = "S_{}_{}_win_{}".format(i,winsize,e)
            #FD for HAT
            data = WD[winname]
            freq,counts = FreqDist(data,bins)
            wFD[winname] = freq
            wFDc[winname] = counts
            #FD for S_HAT
            data = WD[S_winname]
            freq,counts = FreqDist(data,bins_S)
            wFD[S_winname] = freq
            wFDc[S_winname] = counts
                
file = open(analysispath + "WFDFile", "wb")
pickle.dump(wFD,file)
file.close()
                
#Create Folder to store plots of Active Corridor Selection
foldername = "WinFD"
plotpath = analysispath + foldername + "\\"

#Change path for plot storage
if not os.path.exists(plotpath):
    os.makedirs(plotpath)

#Change Directory for storage of analysis files
os.chdir(plotpath)

#USER SET PLOTTING PREFERENCES HERE
Plot = False

#Make plots of Frequency Distributions of Dimensionless Bed Elevation
for i in range(len(reachlist)):
    if Plot == True:
        sectnum = "mn_{}".format(i)
        NumSections = sect[sectnum]
        for w in range(len(Sizes)):
            winsize = Sizes[w]
            number = np.ceil(NumSections/winsize)
            start = int(np.floor(winsize/2))
            end = (NumSections-start)
            extent = np.arange(start,end,1)
            for e in extent:
                winname = "{}_{}_win_{}".format(i,winsize,e)
                S_winname = "S_{}_{}_win_{}".format(i,winsize,e)
                fig, ax = plt.subplots(1, 1)
                freq = wFD[S_winname]
                ax.plot(bins_S, freq, 'r.')
                plt.title(winname)
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel("Dimensionless Bed Elevation")
                plt.ylabel("Frequency")
                plt.ylim(0.001,1)
                plt.savefig(S_winname)
    else:
        break

#USER SET PLOTTING PREFERENCES HERE
Plot = False

#Make plots of Frequency Distribution for Absolute Bed Elevation
for i in range(len(reachlist)):
    if Plot == True:
        sectnum = "mn_{}".format(i)
        NumSections = sect[sectnum]
        for w in range(len(Sizes)):
            winsize = Sizes[w]
            number = np.ceil(NumSections/winsize)
            start = int(np.floor(winsize/2))
            end = (NumSections-start)
            extent = np.arange(start,end,1)
            for e in extent:
                winname = "{}_{}_win_{}".format(i,winsize,e)
                S_winname = "S_{}_{}_win_{}".format(i,winsize,e)
                fig, ax = plt.subplots(1, 1)
                freq = wFD[winname]
                ax.plot(bins, freq, 'r.')
                plt.title(winname)
#                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel("Height Above Thalweg (m)")
                plt.ylabel("Frequency")
                plt.ylim(0.008,0.5)
                plt.savefig(S_winname)
    else:
        break


#%%
#Fit gamma distribution and plot

#Create Folder to store plots of Active Corridor Selection
foldername = "FitGam"
plotpath = analysispath + foldername + "\\"

#Change path for plot storage
if not os.path.exists(plotpath):
    os.makedirs(plotpath)

#Change Directory for storage of analysis files
os.chdir(plotpath)

#USER SET WHICH WINDOW SIZES TO FIT
Sizes = [25]

#USER SET PLOTTING PREFERENCES HERE
Plot = True

for i in range(len(reachlist)):
    if Plot == True:
        sectnum = "mn_{}".format(i)
        NumSections = sect[sectnum]
        for w in range(len(Sizes)):
            winsize = Sizes[w]
            number = np.ceil(NumSections/winsize)
            start = int(np.floor(winsize/2))
            end = (NumSections-start)
            extent = np.arange(start,end,1)
            A_Fits = {}
            for e in extent:
                S_winname = "S_{}_{}_win_{}".format(i,winsize,e)
                freq = wFD[S_winname]
                alpha = np.linspace(0.2,6,num = 101)
                loc = 0
                scale = 1
                a_par,ai,MSE = MinError(freq,bins_S,alpha,loc,scale)
                A_Fits[S_winname] = a_par,MSE
                fig, ax = plt.subplots(1, 1)
                ax.plot(bins_S, freq, 'r.')
                GamDist = gamma.pdf(bins_S,a_par,loc = loc,scale = scale)
                ax.plot(bins_S,GamDist,'b-')
                plt.title(winname)
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel("Height Above Thalweg (m)")
                plt.ylabel("Frequency")
                plt.ylim(0.001,1)
                plt.savefig(S_winname)
    else:
        break