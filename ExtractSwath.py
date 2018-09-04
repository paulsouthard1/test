# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 16:45:49 2018

@author: ps29626
"""

import os
import arcpy
import numpy as np
import pickle
import argparse



def ReadData(lookuppath):
    #Read in pickle file with number of sections
    print("Loading Metadata Pickle File")
    sect = pickle.load(open(lookuppath + "sectnum","rb"))
    #Assign number of sections to variable
    NumSect = sect["NumSections"]
    print(NumSect)
    return NumSect

def LoadLines(lookuppath,storepath,arcpath,demname,NumSect,perform):
    if perform = True:
        #Iterate through number of segments generated
        startx_field = "Topx"
        starty_field = "Topy"
        endx_field = "Botx"
        endy_field = "Boty"
        prj = arcpath + filein + "_fdemasc.prj"
        for j in np.arange(1,NumSect+1):
            #Create key to find text file
            lookup = "CS_{}".format(j) + ".txt"
            print("Creating Line from "+ lookup)
            #Assign parameters for ArcTool
            in_table = lookuppath + "CS_{}".format(j)+".txt"
            out_featureclass = storepath+"CS_{}".format(j)+".shp"
            startx_field = "Topx"
            starty_field = "Topy"
            endx_field = "Botx"
            endy_field = "Boty"
            prj = arcpath + demname + "_fdemasc.prj"
            #Run ArcTool
            arcpy.XYToLine_management(in_table, out_featureclass, startx_field, starty_field, endx_field, endy_field, "GEODESIC", "", prj)
    else:
        print("Lines not loaded")
def SwathPull(lookuppath,storepath,dataset,dataname,NumSect):
    #Read in pickle file with points
    print("Loading points Pickle File")
    points = pickle.load(open(lookuppath + "points","rb"))
    print("Successfully loaded points")
    for j in np.arange(2,NumSect+1):
        key2 = "seg_{}".format(j)
        key1 = "seg_{}".format(j-1)
        print("Extracting swath between "+key1+" and "+key2)
        line1 = points[key1]
        line2 = points[key2]
        maxx = np.max([line1[0],line2[0]])
        minx = np.min([line1[0],line2[0]])
        maxy = np.max([line1[1],line2[1]])
        miny = np.min([line1[2],line2[2]])
        if maxx == minx or maxy == miny:
            print("Invalid swath "+key1)
        else:
            rectangle = str(minx)+" "+str(miny)+" "+str(maxx)+" "+str(maxy)
            outraster = storepath+dataname+"_"+str(j)+".tif"
            arcpy.Clip_management(dataset,rectangle,outraster,"#","-9999","NONE","MAINTAIN_EXTENT")
            out_ascii_file = lookuppath+dataname+"_"+key1+".txt"
            arcpy.RasterToASCII_conversion(outraster,out_ascii_file)
            
def Main(folder,demin,dataset,perform):
    dir_path = os.path.dirname(os.path.realpath(__file__))+"\\"
    arcpath = dir_path + "ArcFiles\\"
    crosspath = dir_path + "CrossFiles\\"
    if not os.path.exists(crosspath):
        os.makedirs(crosspath)
    storepath = arcpath + folder + "\\"
    if not os.path.exists(storepath):
        os.makedirs(storepath)
    lookuppath = crosspath + folder + "\\"
    demname = os.path.basename(demin)
    demname = demname[:demname.find(".")]
#    sr = arcpy.SpatialReference("NAD 1983 UTM Zone 12N")
    NumSect = ReadData(lookuppath)
    LoadLines(lookuppath,storepath,arcpath,demname,NumSect,perform)
    dataname = os.path.basename(dataset)
    dataname = dataname[:dataname.find(".")]
    SwathPull(lookuppath,storepath,dataset,dataname,NumSect)
    
parser = argparse.ArgumentParser(description='Extract swatchs of data along a reach and optionally make shapefiles of segments separating swaths')
parser.add_argument('folder', help='Folder to store resulting files in A.K.A. name describing analysis')
parser.add_argument('demin', help='DEM of Region')
parser.add_argument('dataset', help='dataset to extract swaths from')


args = parser.parse_args()    

Main(folder,demin,dataset,perform)