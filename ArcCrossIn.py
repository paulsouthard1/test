
# coding: utf-8


import os
import arcpy
from arcpy import env
import pickle
import shutil
import argparse


#Putting a new comment






dir_path = os.path.dirname(os.path.realpath(__file__))+"\\"
arcpath = dir_path + "ArcFiles\\"
crosspath = dir_path + "CrossFiles\\"
if not os.path.exists(crosspath):
    os.makedirs(crosspath)
#USER SET FILE NAME HERE - Assign name to this study reach
foldname = "Unveg_1"
storepath = arcpath + foldname + "\\"
if not os.path.exists(storepath):
    os.makedirs(storepath)
lookuppath = crosspath + foldname + "\\"

env.workspace = storepath

#USER SET FILENAME HERE - Define file names
filein = "BE_08"

#Read in reach file to determine # of reaches
ReachFile = lookuppath + "reach_" + foldname + ".txt"
with open(ReachFile,'r') as f:
    #Variable reachlist contains coordinates for each reach to be analyzed
    reachlist = f.readlines()

fullfilename = lookuppath + "sectnum"
#Read in pickle file with number of sections
sect = pickle.load(open(fullfilename,"rb"))

sr = arcpy.SpatialReference("NAD 1983 UTM Zone 12N")
for i in range(len(reachlist)):
    sectnum = "mn_{}".format(i)
    for j in range(sect[sectnum]+1):
        lookup = "CS_{}_{}".format(i,j+1) + ".txt"
        in_table = lookuppath + "CS_{}_{}".format(i,j+1)+".txt"
        out_featureclass = "CS_{}_{}".format(i,j+1)+".shp"
        startx_field = "RBX"
        starty_field = "RBY"
        endx_field = "LBX"
        endy_field = "LBY"
        prj = arcpath + filein + "_fdemasc.prj"
        arcpy.XYToLine_management(in_table, out_featureclass, startx_field, starty_field, endx_field, endy_field, "GEODESIC", "", prj)
        in_line_features = storepath + out_featureclass
        profile_targets = dir_path + filein + "\\" + filein + ".tif"
        out_table = "T_" + lookup
        arcpy.StackProfile_3d(in_line_features, profile_targets, out_table)
        shutil.move(storepath + out_table, lookuppath + out_table)