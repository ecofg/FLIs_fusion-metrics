#! python2
# coding=UTF-8
'''
The total information volume  based on Huffman trees optimal coding.
The total information volume of raster data.
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
inputpath=r" "#Target raster path
outpath=r" "#Output result path
rs=1#resolution
r =30/rs# the redius of the moving window
sds=4#four significant digits after the decimal point

#Environment variable
arcpy.CheckOutExtension("spatial")
arcpy.env.parallProcessingFactor = "50%"
arcpy.gp.overwriteOutput = True
inRas = Raster(inputpath)
desc = arcpy.Describe(inRas)
extent1 = desc.extent
print(extent1)
arcpy.env.extent = extent1
sourceSR = arcpy.Describe(inRas).spatialReference
arcpy.env.outputCoordinateSystem = sourceSR
lowerLeft = arcpy.Point(inRas.extent.XMin, inRas.extent.YMin)
cellSize = inRas.meanCellWidth
print(lowerLeft, cellSize)
gc.set_threshold(100, 5, 5)

weishu=range(1,sds+1)
chengshu=map(lambda x:pow(10,x),weishu)

class Node:
    def __init__(self, freq):
        self.left = None
        self.right = None
        self.father = None
        self.freq = freq

    def __repr__(self):
        return "Node({0.freq!r})".format(self)

    def isLeft(self):
        return self.father.left == self

def creatNodes(freqs):
    return [Node(freq) for freq in freqs]

# creat Huffman-tree
def createHuffmanTree(nodes):
    queue = nodes[:]
    while len(queue) > 1:
        queue.sort(key=lambda item: item.freq)
        node_left = queue.pop(0)
        node_right = queue.pop(0)
        node_father = Node(node_left.freq + node_right.freq)
        node_father.left = node_left
        node_father.right = node_right
        node_left.father = node_father
        node_right.father = node_father
        queue.append(node_father)
    queue[0].father = None
    return queue[0]

# Huffman coding
def huffmanEncoding(nodes, root):
    codes = [""] * len(nodes)
    for i in range(len(nodes)):
        node_tmp = nodes[i]
        while node_tmp != root:
            if node_tmp.isLeft():
                codes[i] = "0" + codes[i]
            else:
                codes[i] = "1" + codes[i]
            node_tmp = node_tmp.father
    return codes
if __name__ == '__main__':
    # normalization
    raster = inRas
    maxValue = raster.maximum
    print"maxvalue:" + str(maxValue)
    minValue = raster.minimum
    print"minvalue:" + str(minValue)
    norraster = (raster - float(minValue)) / (float(maxValue) - float(minValue))
    del raster, gc.garbage[:]
    gc.collect()
    numraster = arcpy.RasterToNumPyArray(norraster)
    dfraster = pd.DataFrame(numraster)
    del norraster, numraster
    arrcomp = np.array(dfraster).flatten().astype(float)
    arrready = arrcomp[-np.isnan(arrcomp)]
    del arrcomp, dfraster, gc.garbage[:]
    gc.collect()
    #  Total information volume
    xinxiliang = 0
    for wei, cheng in zip(weishu, chengshu):
        print("wei is %d" % wei)
        print("cheng is %d" % cheng)
        rastercheng = np.rint(np.around(arrready, decimals=wei) * cheng)
        rasterlist = rastercheng.tolist()
        del rastercheng, gc.garbage[:]
        gc.collect()
        changdu = len(rasterlist)
        print(changdu)
        listcount = Counter(rasterlist)
        del rasterlist, gc.garbage[:]
        gc.collect()
        chars = listcount.keys()
        freqs = listcount.values()
        chars_freqs = zip(chars, freqs)
        nodes = creatNodes([item[1] for item in chars_freqs])
        root = createHuffmanTree(nodes)
        codes = huffmanEncoding(nodes, root)
        chars_codes = zip(chars_freqs, codes)
        chars_codes.sort(key=lambda item: item[0][1])
        sum = 0
        for item in chars_codes:
            sum += item[0][1] * len(item[1])
        xinxiliang += int(math.pow(0.5, wei) / (1 - math.pow(0.5, weishu[-1])) * sum)  # 总信息量,为了加快运算取整。
    print("The total information volume of ** is %f Mbit"%xinxiliang)


