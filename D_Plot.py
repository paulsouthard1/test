# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 14:46:15 2018

@author: ps29626
"""

#%%

import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
from scipy.stats import gamma

direc = os.path.dirname(os.path.realpath(__file__)) + "\\"

#Load Data from pickle files
Elevs = pickle.load(open("efile","rb"))
Dists = pickle.load(open("dfile","rb"))

#Create function to plot frequency distribution of values
def FreqDist(data,bins):
    from numpy import zeros,sum
    counts = np.zeros(len(x))
    for i in range(len(x)):
        if i == 0:
            lower = 0
        else:
            lower= x[i-1]
        upper = x[i]
        for j in range(len(data)):
            if (data[j] > lower) & (data[j] <= upper):
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
            if freq[j] == 0:
                SqErr[j] == 0
            else:
                SqErr[j] = ((GamDist[j]-freq[j])/freq[j])**2
        MSE[i] = mean(SqErr)
    ai = argmin(MSE)
    a_par = alpha[ai]
    return a_par,ai,MSE

#Create loop to produce plots for all 
Windows = [9,15,21]

WDists = {}

#Blank Dictionaries for Raw Height above thalwegs and normalized Height above thalweg
HATs = {}
S_HATs = {}

#Read in starting point for reach
ReachFile = direc + "reach.txt"
with open(ReachFile,'r') as f:
    #Variable reachlist contains coordinates for each reach to be analyzed
    reachlist = f.readlines()

# e profile data X-sect by X-sect
for i in range(len(reachlist)):
    sectnum = "mn_{}".format(i)
    for j in range(49):
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

begin = 26
sectnum = 49
for i in range(len(Windows)):
    winsize = Windows[i]
    number = np.ceil(sectnum/winsize)
    start = int(np.floor(winsize/2))+26
    end = (sectnum-start)+26
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
        WDists[S_winname] = S_Ewin

x = np.logspace(np.log10(0.01),np.log10(30))
alpha = np.linspace(0.4,10,97)

#for i in range(len(Windows)):
#    winsize = Windows[i]
#    number = np.ceil(sectnum/winsize)
#    start = int(np.floor(winsize/2))
#    end = sectnum-start
#    extent = np.arange(start,end,1)
#    for j in extent:
#        S_winname = "S_{}_set_{}".format(winsize,j)
#        data = WDists[S_winname]
#        savename = S_winname+".png"
#        freq,counts = FreqDist(data,x)
#        a_par,ai,MSE = MinError(freq,x,alpha)
#        distline = gamma.pdf(x,a_par,scale = 1.9)
#        fig, ax = plt.subplots(1, 1)
#        ax.plot(x, freq, 'r.')
#        ax.plot(x, distline, 'b-', lw=2, label='best')
#        ax.legend(frameon=False)
#        plt.title(S_winname)
#        plt.xscale('log')
#        plt.yscale('log')
#        plt.xlabel("Dimensionless Bed Elevation")
#        plt.ylabel("Frequency")
#        a_label = "alpha = " + str(a_par)
#        ax.text(0.01, 0.35, a_label, fontsize=15)
#        plt.ylim(0.001,1)
#        plt.savefig(savename)
        
begin = 26
sectnum = 49
for i in range(len(Windows)):
    winsize = Windows[i]
    number = np.ceil(sectnum/winsize)
    start = int(np.floor(winsize/2))+26
    end = (sectnum-start)+26
    extent = np.arange(start,end,1)
    for j in extent:
        S_winname = "S_{}_set_{}".format(winsize,j)
        data = WDists[S_winname]
        savename = S_winname+".png"
        freq,counts = FreqDist(data,x)
        a_par,ai,MSE = MinError(freq,x,alpha)
        distline = gamma.pdf(x,a_par,scale = 1.6)
        fig, ax = plt.subplots(1, 1)
        ax.plot(x, freq, 'r.')
        ax.plot(x, distline, 'b-', lw=2, label='best')
        ax.legend(frameon=False)
        plt.title(S_winname)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel("Dimensionless Bed Elevation")
        plt.ylabel("Frequency")
        a_label = "alpha = " + str(a_par)
        ax.text(0.01, 0.35, a_label, fontsize=15)
        plt.ylim(0.001,1)
        plt.savefig(savename)
        
        
        



