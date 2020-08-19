#! python2
# coding=UTF-8
'''
The fusion PR indicator  based on G1,G2 LULC
Moving window method
'''

import arcpy
from arcpy.sa import *
from arcpy.da import *
import os
import gc
import math
import numpy as np
import pandas as pd
from collections import Counter

#Input parameters
inputg1path=r""#Target raster path, G1 LULC
inputg2path=r""# G2 LULC
outpath=r""#Output result path file(.gdb)

rs=1#resolution
r =30/rs# moving window radius
sds=4#four significant digits after the decimal point

#Environment variable
arcpy.CheckOutExtension("spatial")
arcpy.env.parallProcessingFactor = "50%"
arcpy.gp.overwriteOutput = True
inRas = Raster(inputg1path)
desc = arcpy.Describe(inRas)
extent1 = desc.extent
arcpy.env.extent = extent1
sourceSR = arcpy.Describe(inRas).spatialReference
arcpy.env.outputCoordinateSystem = sourceSR
lowerLeft = arcpy.Point(extent1.XMin, extent1.YMin)
cellSize = inRas.meanCellWidth
print(lowerLeft, cellSize)
gc.set_threshold(100, 5, 5)

# background window
dfno = pd.DataFrame(index=range(2 * r + 1), columns=range(2 * r + 1))
for i in range(2 * r + 1):
    for n in range(2 * r + 1):
        if pow(r - i, 2) + pow(r - n, 2) <= pow(r, 2):  # 圆形
            dfno.iat[i, n] = 1
def className(input):
    numarr=arcpy.RasterToNumPyArray(Raster(input))
    dic=Counter(numarr.flatten().tolist())
    classname=dic.keys()
    del dic,numarr
    return classname
def nuMclass(g1n,g2n):
    list=[]
    for n1 in g1n:
        num=0
        for n2 in g2n:
            if str(n2)[0]==str(n1):
                num+=1
        list.append(num)
    return list

def piWindow(dfclp):
    dfwindow = dfclp.rename(index=dict(map(lambda x, y: [x, y], dfclp.index.tolist(), dfno.index.tolist())),
                            columns=dict(map(lambda x, y: [x, y], dfclp.columns.tolist(), dfno.columns.tolist())))
    dfcomp = dfno * dfwindow
    arrcomp = np.array(dfcomp).flatten().astype(float)
    glist = arrcomp[-np.isnan(arrcomp)].tolist()
    gdic = Counter(glist)
    del arrcomp, dfcomp, dfwindow
    for k, v in gdic.items():
        gdic[k] = round(float(gdic[k])/len(glist),3)
    return gdic
if __name__ == '__main__':
    g1class=className(inputg1path)
    g2class=className(inputg2path)
    numberlist=nuMclass(g1class,g2class)
    print(g2class)
    print(g1class)
    print(numberlist)
    numraster = arcpy.RasterToNumPyArray(inRas)
    dfraster = pd.DataFrame(numraster)
    dflower=pd.DataFrame(arcpy.RasterToNumPyArray(Raster(inputg2path)))
    del numraster,inRas,gc.garbage[:]
    gc.collect()
    hang = dfraster.shape[0]  # 行数
    lie = dfraster.shape[1]  # 列数
    dfzero = pd.DataFrame(np.zeros((hang, lie)))
    # 创建圆形移动窗口
    n = 0
    for l in range(lie - 2 * r):
        print("The progress rate is %3f" % (float(l) / (lie - 2 * r)))
        for h in range(hang - 2 * r):
            n += 1
            # print ("No. is "+ str(n))
            dfg1clp = dfraster.iloc[h:h + 2 * r + 1, l:l + 2 * r + 1]
            dfg2clp = dflower.iloc[h:h + 2 * r + 1, l:l + 2 * r + 1]
            pig1=piWindow(dfg1clp)
            pig2=piWindow(dfg2clp)
            sumg1=0
            nums=0
            for g1,num in zip(g1class,numberlist):
                nums+=num
                sumg2=0
                if pig1[g1]==0:
                    continue
                elif num==1:
                    sumg1+=1
                else:
                    for g2 in g2class[nums-num:nums]:
                        if pig2[g2]!=0:
                            sumg2+=pig2[g2]*math.log(pig2[g2])
                    sumg1+=1-pig1[g1]*sumg2/math.log(num)
            dfzero.iat[h + r, l + r] = sumg1
    arrayzero = np.array(dfzero)
    fusionRaster = arcpy.NumPyArrayToRaster(arrayzero, lowerLeft, cellSize)
    arcpy.DefineProjection_management(fusionRaster, sourceSR)
    fusionRaster.save(outpath)
    print ("FLI_PR is ok")
