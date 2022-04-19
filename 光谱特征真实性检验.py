import numpy as np
import math
from osgeo import gdal


def getMax(x, y):
    maxXs = []
    maxYs = []
    maxXs.append(x[0])
    maxYs.append(y[0])
    for index in range(1, x.__len__() - 1):
        if (y[index] - y[index - 1] > y[index + 1] - y[index]):
            maxXs.append(x[index])
            maxYs.append(y[index])

    maxXs.append(x[x.__len__() - 1])
    maxYs.append(y[y.__len__() - 1])
    return maxXs, maxYs


def getStart(maxXs, maxYs):
    maxIndex = np.argmax(maxYs)
    startPointX = maxXs[maxIndex]
    startPointY = maxYs[maxIndex]
    return maxIndex, startPointX, startPointY


def getNextR(maxYs, maxXs, maxIndex, startPointX, startPointY):
    xNextR = []
    yNextR = []
    pointIndex = maxIndex
    indexPointX = startPointX
    indexPointY = startPointY
    while (pointIndex < maxXs.__len__() - 1):
        kMax = -999999
        for index in range(pointIndex + 1, maxXs.__len__()):
            k = (maxYs[index] - indexPointY) / (maxXs[index] - indexPointX)
            if (kMax < k):
                kMax = k
                nextIndex = index
        pointIndex = nextIndex
        indexPointX = maxXs[pointIndex]
        indexPointY = maxYs[pointIndex]
        xNextR.append(indexPointX)
        yNextR.append(indexPointY)

    if (not xNextR) or ((xNextR[xNextR.__len__() - 1] != maxXs[maxXs.__len__() - 1])):
        xNextR.append(maxXs[maxXs.__len__() - 1])
        yNextR.append(maxYs[maxYs.__len__() - 1])
    return xNextR, yNextR


def getNextL(maxYs, maxXs, maxIndex, startPointX, startPointY):
    xNextL = []
    yNextL = []
    pointIndex = maxIndex
    indexPointX = startPointX
    indexPointY = startPointY
    while (pointIndex > 1):
        kMin = 999999
        for index in range(pointIndex - 1, 0, -1):
            k = (maxYs[index] - indexPointY) / (maxXs[index] - indexPointX)
            if (kMin > k):
                kMin = k
                nextIndex = index
        pointIndex = nextIndex
        indexPointX = maxXs[pointIndex]
        indexPointY = maxYs[pointIndex]
        xNextL.append(indexPointX)
        yNextL.append(indexPointY)

    if (not xNextL) or (xNextL[xNextL.__len__() - 1] != maxXs[0]):
        xNextL.append(maxXs[0])
        yNextL.append(maxYs[0])

    return xNextL, yNextL


def doExtend(xNextL, yNextL, xNextR, yNextR, startPointX, startPointY):
    xNextL.reverse()
    yNextL.reverse()
    xNextL += [startPointX]
    xNextL += xNextR
    yNextL += [startPointY]
    yNextL += yNextR
    return xNextL, yNextL


def doInsert(spcx, xNext, yNext):
    yInsert = []
    xInsert = []
    for i in range(xNext.__len__() - 1):
        x = xNext[i]
        y = yNext[i]
        yInsert.append(y)
        xInsert.append(x)
        k = (yNext[i + 1] - yNext[i]) / (xNext[i + 1] - xNext[i])
        b = yNext[i] - xNext[i] * k
        tup = np.where(spcx == x)
        index = tup[0]
        while x < xNext[i + 1]:
            index += 1
            x = spcx[index]
            if x < xNext[i + 1]:
                y = k * x + b
                xInsert.append(x)
                yInsert.append(y)

    xInsert.append(xNext[xNext.__len__() - 1])
    yInsert.append(yNext[yNext.__len__() - 1])
    return yInsert


def getResult(x, y, yInsert):
    yResult = []
    for i in range(x.__len__()):
        yResult.append(y[i] / yInsert[i])

    return yResult

# def writeResult(flieNmae, x, yResult):
#     fp = open(flieNmae, 'w')
#     for i in range(x.__len__()):
#         fp.write("%.6f  " % x[i])
#         for j in range(9):
#             fp.write("%.6f  " % yResult[j][i])
#         fp.write('\n')


def startDeal(x, y, row, col):
    yResult = np.zeros([row,col])
    yInsert = []
    for i in range(0, col):#####9要改
        maxXs, maxYs = getMax(x, y[:,i])#获取最大值
        maxIndex, startPointX, startPointY = getStart(maxXs, maxYs)#选取极值点
        xNextR, yNextR = getNextR(maxYs, maxXs, maxIndex, startPointX, startPointY)#选取右侧满足要求的极值点
        xNextL, yNextL = getNextL(maxYs, maxXs, maxIndex, startPointX, startPointY)#选取左侧满足要求的极值点
        xNext, yNext = doExtend(xNextL, yNextL, xNextR, yNextR, startPointX, startPointY)#合并选择的点
        Insert = doInsert(x, xNext, yNext)#求得包络线
        yInsert.append(Insert)#包络线
        yResult[:,i] = getResult(x, y[:,i], Insert)#包络线去除

    return y,yInsert,yResult

#光谱相关系数
def speccoff(input,true,row,col):
    #input为待检验光谱，第一列为波长
    #true为参考光谱，第一列为波长
    #row为数据行数，光谱数；col待测组数+1，因为input和true第一列为波长
    SCM = np.zeros(col - 1)#相关系数
    # SCM2 = np.zeros(col - 1)
    for i in range(0,col-1):
        specr = true[:,i+1]  #参考光谱
        spect = input[:,i+1] #验证光谱
        # rmean = np.mean(specr)
        # tmean = np.mean(spect)
        # r = specr - rmean
        # t = spect - tmean
        # x1 = np.sum(r * t)
        # x2 = np.sum(r**2) * np.sum(t**2)
        # SCM[i-1] = x1 / np.sqrt(x2)
        SCM1 = np.corrcoef(specr, spect)
        SCM[i] = SCM1[1][0]
    return SCM

#光谱角
def specsam(input,true,row,col):
    # input为待检验光谱，第一列为波长
    # true为参考光谱，第一列为波长
    # row为数据行数，光谱数；col待测组数+1，因为input和true第一列为波长
    SAM = np.zeros(col - 1)#光谱角匹配
    # SCM2 = np.zeros(col - 1)
    for i in range(0, col - 1):
        specr = true[:, i+1]  # 参考光谱
        spect = input[:, i+1]  # 验证光谱
        x1 = np.sum(specr * spect)
        sumr = np.sum(specr ** 2)
        sumt = np.sum(spect ** 2)
        x2 = np.sqrt(sumr) * np.sqrt(sumt)
        SAM[i] = math.acos(x1/x2)
         #= 1/x3
    return SAM


#连续统去除
def envelope(wn,spcheck, sptrue, row, col):
    #wn为波长数组，spcheck为待检测光谱，sptrue为参考光谱
    #row为数据行数，光谱数；col待测组数
    local = np.zeros(col)#吸收位置
    depth = np.zeros(col)#吸收深度
    coff = np.zeros(col)#相关系数
    sam = np.zeros(col)#光谱角

    # spcheck = input[:,1:col]
    # sptrue = true[:, 1:col]
    [spc, c1, spce] = startDeal(wn, spcheck, row, col)#连续统去除程序
    [spt, t1, spet] = startDeal(wn, sptrue, row, col)
    wn = wn[:, np.newaxis]
    inpute = np.hstack((wn, spce))#为使用speccoff和specsam函数运算合并波长与光谱数组
    truee = np.hstack((wn, spet))

    #连续统去除后计算吸收位置、深度、相关系数、光谱角
    for i in range(col-1):
        depth[i] = np.min(spce[:, i])
        tu = np.where(spce[:, i] == depth[i])
        for st in tu[0]:
            break
        local[i] = wn[int(st)]
        coff[i] = speccoff(inpute[:, [0, i+1]], truee, row, 2)
        sam[i] = specsam(inpute[:, [0, i+1]], truee, row, 2)
    # for i in range(6):
    #     # plt.plot(wn, spc[:,i], color='r')
    #     # plt.plot(wn, c1[i], color='g')
    #     plt.plot(wn, spce[i], color='b')
    #     plt.xlim(0)
    #     plt.ylim(0)
    #     plt.show()
    return  local, depth, coff, sam



def specchatrue(inputfile, truefile, savefile, wstart, wend, shpfile=False, inputwarpfile1=False, inputwarpfile2=False ):
    if not shpfile:
        tshp = 0
    else:
        tshp = 1
    # 读取为np数组
    # input = np.loadtxt(inputfile,  comments=['Column', 'ENVI'])
    # true = np.loadtxt(truefile,  comments=['Column', 'ENVI'])
    # wn = input[:, 0]
    input = gdal.Open(inputfile)
    true = gdal.Open(truefile)
    if tshp == 1:
        inputclass = gdal.Warp(inputwarpfile1, input, format='GTiff', cutlineDSName=shpfile, cropToCutline=True,
                               dstNodata=0)
        trueclass = gdal.Warp(inputwarpfile2, true, format='GTiff', cutlineDSName=shpfile, cropToCutline=True,
                              dstNodata=0)
        array_input = inputclass.ReadAsArray()
        array_true = trueclass.ReadAsArray()
    else:
        array_input = input.ReadAsArray()
        array_true = true.ReadAsArray()

    sample1 = np.size(array_true, 2)
    line1 = np.size(array_true, 1)
    sample2 = np.size(array_input, 2)
    line2 = np.size(array_input, 1)
    if (sample1 == sample2) and (line1 == line2):
        sample = sample1
        line = line1
    else:
        sample = min(sample1, sample2)
        line = min(line1, line2)
        array_true = array_true[:, 0:line, 0:sample]
        array_input = array_input[:, 0:line, 0:sample]

    # 截取波长范围的数组
    wn = np.zeros(true.RasterCount)
    for i in range(1, true.RasterCount + 1):
        wn[i - 1] = float(true.GetRasterBand(i).GetMetadataItem("wavelength"))
    
    wstart = float(wstart)
    wend = float(wend)
    if (wstart <= 0) or (wstart == ''):
        wstart = wn[0]
    if (wend <= wstart) or (wend == ''):
        wend = wn[-1]

    wn[wn > wend] = -1
    inputsbool = wn >= wstart
    wn = wn[inputsbool]

    inputcut = array_input[inputsbool, :, :]
    truecut = array_true[inputsbool, :, :]
    row = np.size(inputcut, 0)  # 光谱
    col = np.size(inputcut, 1)  # 行
    col2 = np.size(inputcut, 2)  # 列
    local = np.zeros((col, col2))
    depth = np.zeros((col, col2))
    coff = np.zeros((col, col2))
    sam = np.zeros((col, col2))

    # input = input[inputsbool, :]
    # true = true[inputsbool, :]
    # row = np.size(input, 0)  # 计算 X 的行数
    # col = np.size(input, 1)  # 计算 X 的列数
    for i in range(col2):
        input1 = inputcut[:, :, i]
        true1 = truecut[:, :, i]
        input1[input1 <= 0] = 1
        true1[true1 <= 0] = 1
        [local[:, i], depth[:, i], coff[:, i], sam[:, i]] = envelope(wn, input1, true1, row, col)
    # local = local[np.newaxis, :]
    # depth = depth[np.newaxis, :]
    # coff = coff[np.newaxis, :]
    # sam = sam[np.newaxis, :]

    # savetxtfile = savefile + '.img'
    driver = gdal.GetDriverByName("ENVI")
    # filename = savefile + '.img'
    outdata = driver.Create((savefile + '.img'), col2, col, 4, gdal.GDT_Float32)
    # outdata.SetGeoTransform(im_geotrans)  # 写入仿射变换参数
    # outdata.SetProjection(im_proj)  # 写入投影
    outdata.GetRasterBand(1).WriteArray(local)
    outdata.GetRasterBand(2).WriteArray(depth)
    outdata.GetRasterBand(3).WriteArray(coff)
    outdata.GetRasterBand(4).WriteArray(sam)

    # with open(savetxtfile, 'ab') as f:
    #     np.savetxt(f, local, fmt='%.6f', delimiter='   ', newline='\n', comments='', header='吸收位置：', encoding='utf-8')
    #     np.savetxt(f, depth, fmt='%.6f', delimiter='   ', newline='\n', comments='', header='吸收深度：', encoding='utf-8')
    #     np.savetxt(f, coff, fmt='%.6f', delimiter='   ', newline='\n', comments='', header='相关系数：', encoding='utf-8')
    #     np.savetxt(f, sam, fmt='%.6f', delimiter='   ', newline='\n', comments='', header='光谱角：', encoding='utf-8')


inputfile = r'C:\Users\agrs\Desktop\to xx\py文件\test\光谱特征真实性检验\gf5_279band'
truefile = r'C:\Users\agrs\Desktop\to xx\py文件\test\光谱特征真实性检验\hymap_to_gf5_279band'
savefile = r'C:\Users\agrs\Desktop\to xx\py文件\test\光谱特征真实性检验\ccc'
# shpfile 默认值=False，如需要，传入地址，对应 app 中的 “检测范围”
shpfile = r'C:\Users\agrs\Desktop\to xx\py文件\test\光谱特征真实性检验\reference.shp'
# inputwarpfile，缓存文件
inputwarpfile1 = r'C:\Users\agrs\Desktop\to xx\example\warp\inputwarp1.tif'
inputwarpfile2 = r'C:\Users\agrs\Desktop\to xx\example\warp\inputwarp2.tif'

# 输入检验波长范围，wstart ~ wend（以nm为单位）
wstart = 800
wend = 1200


specchatrue(inputfile, truefile, savefile, wstart, wend, shpfile, inputwarpfile1, inputwarpfile2)
