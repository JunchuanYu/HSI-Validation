from osgeo import gdal
import numpy as np
import math


# 光谱相关系数
def speccoff(input, true, row, col):
    # input为待检验光谱，第一列为波长
    # true为参考光谱，第一列为波长
    # row为数据行数，光谱数；col待测组数+1，因为input和true第一列为波长
    SCM = np.zeros(col - 1)  # 相关系数
    # SCM2 = np.zeros(col - 1)
    for i in range(0, col - 1):
        specr = true[:, i + 1]  # 参考光谱
        spect = input[:, i + 1]  # 验证光谱
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


# 光谱角
def specsam(input, true, row, col):
    # input为待检验光谱，第一列为波长
    # true为参考光谱，第一列为波长
    # row为数据行数，光谱数；col待测组数+1，因为input和true第一列为波长
    SAM = np.zeros(col - 1)  # 光谱角匹配
    # SCM2 = np.zeros(col - 1)
    for i in range(0, col - 1):
        specr = true[:, i + 1]  # 参考光谱
        spect = input[:, i + 1]  # 验证光谱
        x1 = np.sum(specr * spect)
        sumr = np.sum(specr ** 2)
        sumt = np.sum(spect ** 2)
        x2 = np.sqrt(sumr) * np.sqrt(sumt)
        SAM[i] = math.acos(x1 / x2)
        # = 1/x3
    return SAM


# 欧氏距离
def specdis(input, true, row, col):
    # input为待检验光谱，第一列为波长
    # true为参考光谱，第一列为波长
    # row为数据行数，光谱数；col待测组数+1，因为input和true第一列为波长
    DIS = np.zeros(col - 1)  # 欧氏距离
    # SCM2 = np.zeros(col - 1)
    for i in range(0, col - 1):
        specr = true[:, i + 1]  # 参考光谱
        spect = input[:, i + 1]  # 验证光谱
        t = np.sum((specr - spect) ** 2)
        DIS[i] = np.sqrt(t) / 10000
    return DIS


def reftruecal(inputfile, shpfile, truefile, savefile, inputwarpfile1, inputwarpfile2, wstart, wend, iscoff=False,
               issam=False, isdis=False):
    gdal.AllRegister()

    if not shpfile:
        tshp = 0
    else:
        tshp = 1
    # if not trueshpfile:
    #     if tshp == 1:
    #         trueshpfile = shpfile

    input = gdal.Open(inputfile)
    true = gdal.Open(truefile)

    # shp = shapefile.Reader(shpfile)
    # trueshp = shapefile.Reader(trueshpfile)
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
    # wn1 = true.GetRasterBand(1).GetMetadataItem("wavelength")
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

    # 三个评价指标
    coff = np.zeros([col, col2])
    sam = np.zeros([col, col2])
    dis = np.zeros([col, col2])
    wn = wn[:, np.newaxis]
    for i in range(col2):
        input1 = np.hstack((wn, inputcut[:, :, i]))
        true1 = np.hstack((wn, truecut[:, :, i]))
        if iscoff:
            coff[:, i] = speccoff(input1, true1, row, col + 1)
        if issam:
            sam[:, i] = specsam(input1, true1, row, col + 1)
        if isdis:
            dis[:, i] = specdis(input1, true1, row, col + 1)
    # output = np.zeros(col, col2, 3)
    # coff = coff[:, np.newaxis]
    # sam = sam[:, np.newaxis]
    # dis = dis[:, np.newaxis]
    # output[:, :, 0] = coff
    # output[:, :, 1] = sam
    # output[:, :, 2] = dis
    # a = coff.dtype.name
    driver = gdal.GetDriverByName("ENVI")
    # filename = savefile + '.img'
    outdata = driver.Create((savefile + '.img'), col2, col, 3, gdal.GDT_Float32)

    outdata.GetRasterBand(1).WriteArray(coff)
    outdata.GetRasterBand(2).WriteArray(sam)
    outdata.GetRasterBand(3).WriteArray(dis)


# inputfile 输入待真实性检验文件
inputfile = r'C:\Users\agrs\Desktop\to xx\py文件\test\反射率产品真实性检验\间接\gf5_279band'
# shpfile 默认值=False，如需要，传入地址，对应 app 中的 “检测范围”
shpfile = r'C:\Users\agrs\Desktop\to xx\py文件\test\反射率产品真实性检验\间接\reference.shp'
# truefile 输入参考文件
truefile = r'C:\Users\agrs\Desktop\to xx\py文件\test\反射率产品真实性检验\间接\hymap_to_gf5_279band'
savefile = r'C:\Users\agrs\Desktop\to xx\py文件\test\反射率产品真实性检验\间接\ccc'
# inputwarpfile，默认值=False， 无用的缓存文件，待优化……………………
inputwarpfile1 = r'C:\Users\agrs\Desktop\to xx\example\warp\inputwarp.tif'
inputwarpfile2 = r'C:\Users\agrs\Desktop\to xx\example\warp\inputwarp2.tif'

# 输入检验波长范围，wstart ~ wend（以nm为单位）
wstart = 800
wend = 1000

# 三个评价指标，默认值=Fales（光谱相关系数，光谱角，欧式距离）
iscoff = True
issam = True
isdis = False

reftruecal(inputfile, shpfile, truefile, savefile, inputwarpfile1, inputwarpfile2, wstart, wend, iscoff, issam, isdis)
