#! python2
# coding=UTF-8
'''
The fusion indicators of NP,PD,PR,PRD,SHDI based on Fragstats result.
First, calculate the corresponding indicators of G1,G2 LULC through Fragstats by moving window method.
Second, some data needs to be preprocessed according to detailed requirements.
'''

import arcpy
from arcpy.sa import *
from arcpy.da import *
import os
import gc
import math

#Input parameters
outputfile=r""# output geodatabase
clipFeatures =r""# area boundary
# The landscape path classes level data of NP is stored in geodatabase, named by default.
# Use the pland indicator as area ratio  p(), which still stored in geodatabase.
# Use G1,G2 as prefixes for distinguishing source data respectively.
# Convert the nodata or null value of all source data to 0.
arcpy.env.workspace=r""#stored geodatabase
temporarypath=r""#temporary geodatabase to reduce memory overflow

# classification of LULC, for example:
g1class = ["1","2","3","4","5","6"]
g2class = ['11','12','21','22','23','24','31','32','33','41','42','43','44','45','46','51','52','53','61','62','63','64','65','66','67']
mquantity=[2,4,3,6,3,7]#The max-number of g2class inside corresponding firt grade LULC
A=100 #The area of moving window,according to the actual situation.

#Environment variable
arcpy.CheckOutExtension("spatial")
arcpy.env.parallProcessingFactor = "50%"
arcpy.gp.overwriteOutput = True
gc.set_threshold(100, 10, 10)
arcpy.CreateFileGDB_management(temporarypath,"temporary")
linshi=os.path.join(temporarypath,"temporary")
desc = arcpy.Describe(clipFeatures)
extent1=desc.extent
arcpy.env.extent=extent1
sourceSR = arcpy.Describe(clipFeatures).spatialReference
arcpy.env.outputCoordinateSystem = sourceSR

def fusionNP(linshi):
    sum=[]
    for g1 in g1class:
        npg1 = Raster("G1" + "NP_" + g1)
        weight=Raster(os.path.join(linshi, "wt" + g1))
        npg1part = npg1 * (1 + weight)
        npg1partpath = os.path.join(linshi, "G1g1class" + "NP" + g1)
        npg1part.save(npg1partpath)
        del npg1part,npg1
        sum.append(npg1partpath)
    funp=CellStatistics(sum, "SUM", "DATA")
    arcpy.Clip_management(funp, str(extent1), os.path.join(outputfile, "FLI_NP"), clipFeatures, "",
                          "ClippingGeometry",
                          "NO_MAINTAIN_EXTENT")
    del funp,gc.garbage[:]
    gc.collect()

def fusionPR(linshi):
    sum=[]
    prraster = Raster("G1" + "PR")
    for g1 in g1class:
        weight=Raster(os.path.join(linshi, "wt" + g1))
        prg1=Con(prraster==int(g1),1,0)
        prg1part = prg1 * (1 + weight)
        prg1partpath = os.path.join(linshi, "G1g1class" + "PR" + g1)
        prg1part.save(prg1partpath)
        del prg1part
        sum.append(prg1partpath)
    fupr=CellStatistics(sum, "SUM", "DATA")
    arcpy.Clip_management(fupr, str(extent1), os.path.join(outputfile, "FLI_PR"), clipFeatures, "",
                          "ClippingGeometry",
                          "NO_MAINTAIN_EXTENT")
    del fupr,prraster,gc.garbage[:]
    gc.collect()

def fusionSHDI(linshi):
    sum=[]
    for g1 in g1class:
        g1pland = Raster("G1" + "pland" + "_" + g1)
        weight=Raster(os.path.join(linshi, "wt" + g1))
        shdig1part = -(g1pland / 100) * Ln(g1pland / 100) * (1 + weight)
        shdig1partpath = os.path.join(linshi, "G1g1class" + "SHDI" + g1)
        shdig1part.save(shdig1partpath)
        del shdig1part,g1pland
        sum.append(shdig1partpath)
    fushdi=CellStatistics(sum, "SUM", "DATA")
    arcpy.Clip_management(fushdi, str(extent1), os.path.join(outputfile, "FLI_SHDI"), clipFeatures, "",
                          "ClippingGeometry",
                          "NO_MAINTAIN_EXTENT")
    del fushdi,gc.garbage[:]
    gc.collect()

def fusionPD(linshi):
    sum=[]
    for g1 in g1class:
        pdg1 = Raster("G1" + "PD_" + g1)
        weight=Raster(os.path.join(linshi, "wt" + g1))
        pdg1part = pdg1 * (1 + weight)/A
        pdg1partpath = os.path.join(linshi, "G1g1class" + "PD" + g1)
        pdg1part.save(pdg1partpath)
        del pdg1part,pdg1
        sum.append(pdg1partpath)
    fupd=CellStatistics(sum, "SUM", "DATA")
    arcpy.Clip_management(fupd, str(extent1), os.path.join(outputfile, "FLI_PD"), clipFeatures, "",
                          "ClippingGeometry",
                          "NO_MAINTAIN_EXTENT")
    del fupd,gc.garbage[:]
    gc.collect()

def fusionPRD(linshi):
    sum=[]
    prraster = Raster("G1" + "PRD")
    for g1 in g1class:
        weight=Raster(os.path.join(linshi, "wt" + g1))
        prdg1=Con(prraster==int(g1),1,0)
        prdg1part = prdg1 * (1 + weight)/A
        prdg1partpath = os.path.join(linshi, "G1g1class" + "PRD" + g1)
        prdg1part.save(prdg1partpath)
        del prdg1part
        sum.append(prdg1partpath)
    fuprd=CellStatistics(sum, "SUM", "DATA")
    arcpy.Clip_management(fuprd, str(extent1), os.path.join(outputfile, "FLI_PRD"), clipFeatures, "",
                          "ClippingGeometry",
                          "NO_MAINTAIN_EXTENT")
    del fuprd,prraster,gc.garbage[:]
    gc.collect()

def cleanLinshi(linshi):
    for data in arcpy.ListRasters(linshi):
        arcpy.Delete_management(data)
    gc.collect()

for g1, mq in zip(g1class, mquantity):
    g1pland =Raster("G1" + "pland" + "_" + g1)
    g2tailsum=[]#计算sumg1tail 数据(Atail)
    for g2 in g2class:
        if g2[0] == g1:
            print(g1, g2, g2[0])
            rasterg2 = Raster("G2" + "pland" + "_" + g2)
            pji = rasterg2 / g1pland
            del rasterg2
            gc.collect()
            tail = (-1) * pji * Ln(pji)
            contail=Con(IsNull(tail), 0, tail)
            contailpath=os.path.join(linshi, "G2" + "pland" + "_" + g2 + "contail")
            contail.save(contailpath)
            g2tailsum.append(contailpath)
            del tail,contail,pji,gc.garbage[:]
            gc.collect()
    sumg1tail = CellStatistics(g2tailsum, "SUM", "DATA")
    weightg1=(g1pland / 100) * sumg1tail / math.log(mq)
    weightg1path = os.path.join(linshi, "wt" + g1)
    sumg1tail.save(weightg1path)
if __name__ == '__main__':
    fusionNP(linshi)
    fusionPD(linshi)
    fusionPR(linshi)
    fusionPRD(linshi)
    fusionSHDI(linshi)
    cleanLinshi(linshi)


