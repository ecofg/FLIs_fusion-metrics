#! python2
# coding=UTF-8
'''
The information volume  based on Huffman trees optimal coding.
Moving window method.
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
r =30/rs# moving window radius
sds=4#four significant digits after the decimal point

#Environment variable
arcpy.CheckOutExtension("spatial")
arcpy.env.parallProcessingFactor = "50%"
arcpy.gp.overwriteOutput = True
inRas = Raster(inputpath)  # 提取信息，作为背景信息。
desc = arcpy.Describe(inRas)  # 描述函数
extent1 = desc.extent  # 使用提取矩形的四个顶点坐标
print(extent1)
arcpy.env.extent = extent1  # 环境设置范围，保证所有的栅格处理范围一致，不会发生偏移。
sourceSR = arcpy.Describe(inRas).spatialReference  # 获取原坐标系
arcpy.env.outputCoordinateSystem = sourceSR  # 默认输出坐标系
lowerLeft = arcpy.Point(inRas.extent.XMin, inRas.extent.YMin)  # 左下顶点坐标，为了定位使用。
cellSize = inRas.meanCellWidth
print(lowerLeft, cellSize)
gc.set_threshold(100, 5, 5)

weishu=range(1,sds+1)
chengshu=map(lambda x:pow(10,x),weishu)#需要验证的小数位 10000对应小数后4位

# 创建全为空值的矩形数据帧, background window
dfno = pd.DataFrame(index=range(2 * r + 1), columns=range(2 * r + 1))
# 在建立相应的圆形数据帧
for i in range(2 * r + 1):
    for n in range(2 * r + 1):
        if pow(r - i, 2) + pow(r - n, 2) <= pow(r, 2):  # 圆形
            dfno.iat[i, n] = 1
class Node:
    def __init__(self, freq):
        self.left = None
        self.right = None
        self.father = None
        self.freq = freq

    def __repr__(self):  # 显示属性#change the string representation of instances
        return "Node({0.freq!r})".format(self)

    def isLeft(self):
        return self.father.left == self


# 创建叶子节点  create nodes
def creatNodes(freqs):
    return [Node(freq) for freq in freqs]


# creat Huffman-tree 创建Huffmans树
def createHuffmanTree(nodes):
    queue = nodes[:]  # copy of the nodes
    while len(queue) > 1:
        queue.sort(key=lambda item: item.freq)  # sort the objects by certain attribute
        node_left = queue.pop(0)  # 移除第一个元素
        node_right = queue.pop(0)
        node_father = Node(node_left.freq + node_right.freq)
        node_father.left = node_left
        node_father.right = node_right
        node_left.father = node_father
        node_right.father = node_father
        queue.append(node_father)
    queue[0].father = None  # self.father = None ,可省略
    return queue[0]  # root node


# Huffman 编码
def huffmanEncoding(nodes, root):
    codes = [""] * len(nodes)  # 叶子节点数即为编码数
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
    # 归一化处理
    raster = inRas
    maxValue = raster.maximum
    print"maxvalue:" + str(maxValue)
    minValue = raster.minimum
    print"minvalue:" + str(minValue)
    norraster = (raster - float(minValue)) / (float(maxValue) - float(minValue))
    del raster, gc.garbage[:]
    gc.collect()
    # 转为pd
    numraster = arcpy.RasterToNumPyArray(norraster)
    dfraster = pd.DataFrame(numraster)
    del norraster, numraster
    hang = dfraster.shape[0]  # 行数
    lie = dfraster.shape[1]  # 列数
    dfzero = pd.DataFrame(np.zeros((hang, lie)))  # 创建所有值为0的数据帧,作为填充数据背景
    # 创建圆形移动窗口
    n = 0
    for l in range(lie - 2 * r):
        print("The progress rate is %3f"%(l/(lie-2*r)))
        for h in range(hang - 2 * r):
            n += 1  # 计数起点，矩形框的左上角
            # print ("No. is "+ str(n))
            dfclp = dfraster.iloc[h:h + 2 * r + 1,l:l + 2 * r + 1]
            # 更改行列号，重置，使其可以运算
            dfwindow = dfclp.rename(index=dict(map(lambda x, y: [x, y], dfclp.index.tolist(), dfno.index.tolist())),columns=dict(map(lambda x, y: [x, y], dfclp.columns.tolist(), dfno.columns.tolist())))
            dfcomp = dfno * dfwindow
            arrcomp = np.array(dfcomp).flatten().astype(float)  # 先转一维
            arrready = arrcomp[-np.isnan(arrcomp)]  # 删除空值
            # print(arrready)
            del arrcomp, dfcomp, dfwindow, dfclp
            # 准备进入计算过程，计算信息量，0.5系数法
            xinxiliang = 0
            for wei, cheng in zip(weishu, chengshu):  # 相应的小数位处理
                # print("wei is %d" % wei)
                # print("cheng is %d" % cheng)
                rastercheng = np.rint(np.around(arrready, decimals=wei) * cheng)  # 转为对应的整值
                rasterlist = rastercheng.tolist()
                del rastercheng, gc.garbage[:]
                gc.collect()
                changdu = len(rasterlist)
                # print(changdu)
                listcount = Counter(rasterlist)  # 得到计数的字典
                del rasterlist, gc.garbage[:]
                gc.collect()
                # 最优编码信息量统计
                chars = listcount.keys()
                freqs = listcount.values()
                chars_freqs = zip(chars, freqs)
                nodes = creatNodes([item[1] for item in chars_freqs])
                root = createHuffmanTree(nodes)
                codes = huffmanEncoding(nodes, root)
                chars_codes = zip(chars_freqs, codes)
                chars_codes.sort(key=lambda item: item[0][1])  # sort the result by the freqs
                sum = 0
                for item in chars_codes:
                    # print('Character:%s freq:%-2d encoding:%s' % (item[0][0], item[0][1], item[1]))  # 得到的编码为字符串，可以转数字（十进制），再根据（位数-1）*2+..等转为二进制。
                    # 直接len(item[1])就是其对应的二进制编码的长度，也是其占用的二进制字长，一个二进制位包含的信息为一个比特。所以其长度就是比特数。
                    sum += item[0][1] * len(item[1])  # 只要总结果
                # 权重二分法，分母取其之和。
                xinxiliang += int(math.pow(0.5, wei) / (1 - math.pow(0.5, weishu[-1])) * sum)
            # print("xinxiliang is %f"%xinxiliang)
            # 中新赋值 为 （d+r,d+r),，没有+1，注意语法
            dfzero.iat[h + r, l + r] = xinxiliang
            del arrready, gc.garbage[:]
    arrayzero = np.array(dfzero)
    del dfzero, dfraster, gc.garbage[:]
    newRaster = arcpy.NumPyArrayToRaster(arrayzero, lowerLeft, cellSize)
    print("newRaster is ok")
    # 定义投影则是原来坐标系无或者unknown
    arcpy.DefineProjection_management(newRaster, sourceSR)  # 无坐标系需要定义
    newRaster.save(outpath)
    print ("saving is ok")
    del newRaster, arrayzero, gc.garbage[:]
    gc.collect()

