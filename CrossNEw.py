# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 12:58:02 2018

@author: ps29626
"""

import numpy as np
import os
import pickle
import shutil
import argparse

parser = argparse.ArgumentParser(description='Produce some vertical lines along a flowpath to sample channel characteristics longitudinally')
parser.add_argument('folder', help='Folder to store resulting files in A.K.A. name describing analysis')
parser.add_argument('demin', help='DEM of Region')
parser.add_argument('reachin', help='X,Y Coordinates of reach')
parser.add_argument('--spacing',  type=int, nargs='?', default=50, help='Spacing between lines, length of sampling window')
parser.add_argument('--width',  type=int, nargs='?', default=300, help='Spacing between lines, length of sampling window')
args = parser.parse_args()






#Define neighbor function
def neighbors(im, i, j, d=1):
    b = im[i-d:i+d+1, j-d:j+d+1].flatten()
    # remove the element (i,j)
    n = np.hstack((b[:len(b)//2],b[len(b)//2+1:] ))
    return n

#Define loadascii function
def loadascii(DEMA, FAA):
    print("Loading ASCII file " + DEMA)
    #Load DEM
    DArray=np.loadtxt(DEMA,skiprows=6)
    #Load FAA
    FAArray=np.loadtxt(FAA,skiprows=6)
    return DArray, FAArray

#Define read header function
def readheader(DEMA):
    with open(DEMA,'r') as f:
        print("Reading pertinent header info from " + DEMA)
        #Read in number of columns
        ncols = str(f.readlines(1))
        ncols = ''.join(filter(str.isdigit, ncols))
        #Read in number of rows
        nrows = str(f.readlines(2))
        nrows = ''.join(filter(str.isdigit, nrows))
        #Read in X coordinate of lower left corner
        xllcorner = str(f.readlines(3))
        xllcorner = ''.join(filter(str.isdigit, xllcorner))
        xll = float(xllcorner[:6]+"."+xllcorner[6:])
        #Read in y coordinate of lower left corner
        yllcorner = str(f.readlines(4))
        yllcorner = ''.join(filter(str.isdigit, yllcorner))
        yll = float(yllcorner[:7] + "." + yllcorner[7:])
        #Read in cellsize
        cellsize = str(f.readlines(5))
        cellsize = int(''.join(filter(str.isdigit, cellsize)))
        #Read in cell value for NoData cells
        NODATA = str(f.readlines(6))
        NODATA = ''.join(filter(str.isdigit, NODATA))
        #Find remainders
        xllr = xll%1
        yllr = yll%1
        #Find y coordinate lower left corner of upper left cell for indexing cells
        yul = yll + (float(nrows)-1)
    return ncols, nrows, xll, yll, cellsize, NODATA, xllr, yllr, yul


#Define read reach coordinates function
def readreach(reachin, storepath, folder):
    print("Reading reach beginning and terminus from "+reachin)
    #Read in starting point for reach
    with open(reachin,'r') as f:
        #Variable reachlist contains coordinates for each reach to be analyzed
        reachlist = f.readlines()
        reachlist = reachlist[0]
        shutil.copyfile(reachin,storepath+"reach_"+folder+".txt")
    return reachlist

def pointlist(reachlist,FAArray,xll,yll,yul,xllr,yllr,cellsize,spacing,width,storepath):
    # Create blank dictionary to store profile points along flow accumulation line
    points = {}
    # Create blank dictionary to store profile indices along flow accumulation line
    inds = {}
    #Define exact coordinates of start and end cell corners
    x1 = float(reachlist[0:6])+xllr
    y1 = float(reachlist[7:14])+yllr
    x2 = float(reachlist[15:21])+xllr
    y2 = float(reachlist[22:29])+yllr
    print("Creating cross-sections for reach ("+str(x1)+", "+str(y1)+") to ("+str(x2)+", "+str(y2)+")")
    #Find starting indices
    l1 = round((x1-xll)/cellsize)
    k1 = round((yul-y1)/cellsize)
    #Define initial k,l
    l = l1
    k = k1
    #Find ending indices
    l2 = round((x2-xll)/cellsize)
    k2 = round((yul-y2)/cellsize)
    #Create zero-value variable to test whether spacing distance has been achieved on each iteration
    longthresh = 0
    #Create zero-value variable to order profile points
    seg = 0
    # Create blank variables to store k and l coordinates for each reach
    klist = np.zeros(1000000)
    llist = np.zeros(1000000)
    #Run loop that works through reach using flow accumulation array and selects cell with maximum accumulation
    for m in range(0,999999):
        #Store k's and l's in variable
        klist[m] = k
        llist[m] = l
        #Break loop if at end of reach
        if (k == k2) & (l == l2):
            print("Termination point found")
            break
        #Create a new cross-section point if spacing has been exceeded
        elif longthresh >= spacing:
            #Define point number
            seg = seg + 1
            #Compute and store coordinates of center of cell that is 
            tempx = xll + (l * cellsize) + (cellsize/2)
            tempy = yul - (k * cellsize) + (cellsize/2)
            #Store indices where points were selected
            inds["seg_{}".format(seg)] = [m]
            #Store coordinates to base cross-section line
            topy = tempy + (spacing/2)
            boty = tempy - (spacing/2)
            points["seg_{}".format(seg)] = [tempx,tempy]
            # Save endpoints to a textfile
            savename = "CS_{}".format(seg)+".txt"
            f = open(storepath+savename, 'w')
            f.write("Topx, Topy, Botx, Boty" + "\n")
            f.write(str(tempx) + ", " + str(topy) + ", " + str(tempx) + ", " + str(boty))
            f.close()
#            shutil.move(dir_path+savename,storepath+savename)
            longthresh = 0
            #Run Neighbor function again to get last point
            Temp = neighbors(FAArray,k,l, d=1)
            j = np.argmax(Temp)
            if j == 0:
                k = k-1
                l = l-1
                Long = np.sqrt(2)*cellsize
            elif j == 1:
                k = k-1
                l = l
                Long = cellsize
            elif j == 2:
                k = k-1
                l = l+1
                Long = np.sqrt(2)*cellsize
            elif j == 3:
                k = k
                l = l-1
                Long = cellsize
            elif j == 4:
                k = k
                l = l+1
                Long = cellsize
            elif j == 5:
                k = k+1
                l = l-1
                Long = np.sqrt(2)*cellsize
            elif j == 6:
                k = k+1
                l = l
                Long = cellsize
            else:
                k = k+1
                l = l+1
                Long = np.sqrt(2)*cellsize
            longthresh = longthresh + Long
        else:
            Temp = neighbors(FAArray,k,l, d=1)
            j = np.argmax(Temp)
            if j == 0:
                k = k-1
                l = l-1
                Long = np.sqrt(2)*cellsize
            elif j == 1:
                k = k-1
                l = l
                Long = cellsize
            elif j == 2:
                k = k-1
                l = l+1
                Long = np.sqrt(2)*cellsize
            elif j == 3:
                k = k
                l = l-1
                Long = cellsize
            elif j == 4:
                k = k
                l = l+1
                Long = cellsize
            elif j == 5:
                k = k+1
                l = l-1
                Long = np.sqrt(2)*cellsize
            elif j == 6:
                k = k+1
                l = l
                Long = cellsize
            else:
                k = k+1
                l = l+1
                Long = np.sqrt(2)*cellsize
            # Keep track of how far away the last cross-section is
            longthresh = longthresh + Long
    return m
            
def storemeta(dir_path, storepath, folder, demin, spacing, width, m):
    #Create blank processnotes dictionary
    sect = {}
    #Store number of files for export as pickle file
    sect["NumSections"] = m
    #Save spacing in pickle file
    sect["Spacing"] = spacing
    sect["Width"] = width
    sect["Folder"] = folder
    sect["Reach"] = "reach_"+folder+".txt"
    file = open(dir_path + "sectnum", "wb")
    pickle.dump(sect,file,protocol=2)
    file.close()
    shutil.move(dir_path+"sectnum",storepath+"sectnum")
    
def Main(folder, demin, reachin, spacing, width):
    #Define root directory
    dir_path = os.path.dirname(os.path.realpath(__file__))+"\\"
    #Define Arc File Directory
    arcpath = dir_path + "ArcFiles\\"
    #Define Cross File Directory
    crosspath = dir_path + "CrossFiles\\"
    if not os.path.exists(crosspath):
        os.makedirs(crosspath)
    #Define directory for results storage
    storepath = crosspath + folder + "\\"
    if not os.path.exists(storepath):
        os.makedirs(storepath)
    #Pull out filename
    filein = os.path.basename(demin)
    filein = filein[:filein.find(".")]
    #Define names for ASCII files based on DEM file name
    DEMA = arcpath + filein + "_demasc.txt"
    FAA = arcpath + filein + "_faasc.txt"
    #Load in ASCII's as Arrays
    DArray, FAArray = loadascii(DEMA, FAA)
    #Load in Header information from DEM ASCII
    ncols, nrows, xll, yll, cellsize, NODATA, xllr, yllr, yul = readheader(DEMA)
    #Load in reach endpoints for analysis
    reachlist = readreach(reachin, storepath, folder)
    m = pointlist(reachlist,FAArray,xll,yll,yul,xllr,yllr,cellsize,spacing,width,storepath)
    storemeta(dir_path, storepath, folder, demin, spacing, width, m)
    
#Main(args.folder, args.demin, args.reachin, args.spacing, args.width)
Main(folder, demin, reachin, spacing, width)