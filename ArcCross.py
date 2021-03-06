#%%
# coding: utf-8

import numpy as np
import os
import pickle
import shutil

#Define neighbor function
def neighbors(im, i, j, d=1):
    b = im[i-d:i+d+1, j-d:j+d+1].flatten()
    # remove the element (i,j)
    n = np.hstack((b[:len(b)//2],b[len(b)//2+1:] ))
    return n

#Define mean angles function
from math import sin,cos,atan2,pi
import numpy as np
def meanangle(angles,weights=0,setting='degrees'):
    '''computes the mean angle'''
    if weights==0:
         weights=np.ones(len(angles))
    sumsin=0
    sumcos=0
    if setting=='degrees':
        angles=np.array(angles)*pi/180
    for i in range(len(angles)):
        sumsin+=weights[i]/sum(weights)*sin(angles[i])
        sumcos+=weights[i]/sum(weights)*cos(angles[i])
    average=atan2(sumsin,sumcos)
    if setting=='degrees':
        average=average*180/pi
    return average


dir_path = os.path.dirname(os.path.realpath(__file__))+"\\"
arcpath = dir_path + "ArcFiles\\"
crosspath = dir_path + "CrossFiles\\"
if not os.path.exists(crosspath):
    os.makedirs(crosspath)
#USER SET FILE NAME HERE - Assign name to this study reach
foldname = args.foldname
storepath = crosspath + foldname + "\\"
if not os.path.exists(storepath):
    os.makedirs(storepath)

#%%
#USER SET FILENAME HERE - Define file names
filein = "11_FBE"
DEMA = arcpath + filein + "_demasc.txt"
FAA = arcpath + filein + "_faasc.txt"
FDEMA = arcpath + filein + "_fdemasc.txt"
#USER SET SPACING HERE - Define spacing between profiles
spacing = 20
#USER SET SECTION-LENGTH HERE - Define length of cross section
section = 300

#Load DEM
DArray=np.loadtxt(DEMA,skiprows=6)
#Load FAA
FAArray=np.loadtxt(FAA,skiprows=6)
#Load Filled DEM
FDArray=np.loadtxt(FDEMA,skiprows=6)

#Read in pertinent header info
with open(DEMA,'r') as f:
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

#Read in starting point for reach
ReachFile = dir_path + "reach.txt"
with open(ReachFile,'r') as f:
    #Variable reachlist contains coordinates for each reach to be analyzed
    reachlist = f.readlines()
    
shutil.move(ReachFile,storepath+"reach_"+foldname+".txt")
#%%
# Create blank dictionary to store figure titles
title = {}
# Create blank dictionary to store profile points along flow accumulation line
points = {}
# Create blank dictionary to store profile indices along flow accumulation line
inds = {}
# Create a blank dictionary to store list of k coordinates along flow accumulation line
kcoords = {}
# Create a blank dictionary to store list of l coordinates along flow accumulation line
lcoords = {}
# Create a blank dictionary to store number of sections for each reach to pickle later on
sect = {}

#Define points along reach for analysis, store them as a text file, create shapefile of points
for i in range(len(reachlist)):
    tempreach = reachlist[i]
    # Define exact coordinates of start and end cell corners, and then create title for plots
    x1 = float(tempreach[0:6])+xllr
    y1 = float(tempreach[7:14])+yllr
    x2 = float(tempreach[15:21])+xllr
    y2 = float(tempreach[22:29])+yllr
    title[i] = "Cross section: Reach " + str(i+1)
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
            points["seg_{}".format(seg)] = [tempx,tempy]
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
    #Create list of XY's to create cross-section lines
    for m in range(len(points)):
        # Call index along flow accumulation line where a given cross-section was stored
        seg = m + 1
        Step = inds["seg_{}".format(seg)]
        # 1st of 3 averaged lines
        # Find position of cell directly before cross-section point
        step = Step[0]-1
        backside = [llist[step],klist[step]]
        # Find position of cell directly after cross-section point
        step = step + 2
        frontside = [llist[step],klist[step]]
        #2nd of 3 averaged lines
        # Find position of cell two before cross-section point
        step = Step[0]-5
        backside2 = [llist[step],klist[step]]
        # Find position of cell two after cross-section point
        step = step + 10
        frontside2 = [llist[step],klist[step]]
        # 3rd of 3 averaged lines
        # Find position of cell two before cross-section point
        step = Step[0]-10
        backside3 = [llist[step],klist[step]]
        # Find position of cell two after cross-section point
        step = step + 20
        frontside3 = [llist[step],klist[step]]
        # Average 3 lines on either side to get angle
        #1st theta calc
        # Calculate angle of flow direction from two adjacent points
        if frontside[0] > backside[0]:
            theta = np.arctan((-(frontside[1]-backside[1]))/(frontside[0]-backside[0]))
        else:
            theta = np.arctan((-(frontside[1]-backside[1]))/(-(frontside[0]-backside[0])))
            theta = np.pi - theta
        print(np.degrees(theta))
        #2nd theta calc
        # Calculate angle of flow direction from two adjacent points
        if frontside2[0] > backside2[0]:
            theta2 = np.arctan((-(frontside2[1]-backside2[1]))/(frontside2[0]-backside2[0]))
        else:
            theta2 = np.arctan((-(frontside2[1]-backside2[1]))/(-(frontside2[0]-backside2[0])))
            theta2 = np.pi - theta2
        print(np.degrees(theta2))
        #3rd theta calc
        # Calculate angle of flow direction from two adjacent points
        if frontside3[0] > backside3[0]:
            theta3 = np.arctan((-(frontside3[1]-backside3[1]))/(frontside3[0]-backside3[0]))
        else:
            theta3 = np.arctan((-(frontside3[1]-backside3[1]))/(-(frontside3[0]-backside3[0])))
            theta3 = np.pi - theta3
        print(np.degrees(theta3))
        #Perform averaging
        angles = [theta,theta2,theta3]
        theta = meanangle(angles,weights=0,setting='degrees')
        print(np.degrees(theta))
        print("NEW")
        #Find angle of line perpendicular to stream
        theta = theta - (np.pi/2)
        temppoints = points["seg_{}".format(seg)]
        #Calculate endpoints of 10 m line centered on the flow accumulation line
        RBX = temppoints[0]+(np.cos(theta)*(section/2))
        RBY = temppoints[1]+(np.sin(theta)*(section/2))
        LBX = temppoints[0]-(np.cos(theta)*(section/2))
        LBY = temppoints[1]-(np.sin(theta)*(section/2))
        # Save endpoints to a textfile
        savename = "CS_{}_{}".format(i,seg)+".txt"
        f = open(savename, 'w')
        f.write("RBX, RBY, LBX, LBY" + "\n")
        f.write(str(RBX) + ", " + str(RBY) + ", " + str(LBX) + ", " + str(LBY))
        f.close()
        shutil.move(dir_path+savename,storepath+savename)
    # Store number of files for export as pickle file
    mnum = "mn_{}".format(i)
    sect[mnum] = m
# Save spacing in pickle file
sect["spacing"] = spacing
file = open(dir_path + "sectnum", "wb")
pickle.dump(sect,file,protocol=2)
file.close()
shutil.move(dir_path+"sectnum",storepath+"sectnum")