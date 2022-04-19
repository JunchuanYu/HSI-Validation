import numpy as np
import math
from osgeo import gdal


def calsrf(wn, FWH):
    fwhm = FWH[:, 1]
    wn2 = FWH[:, 0]
    n1 = np.size(wn, 0)
    n2 = np.size(wn2, 0)
    res = np.zeros([n2, n1])
    t = 4
    for i in range(n2):
        sigma = fwhm[i] / math.sqrt(math.log(t))
        sigma = (sigma * sigma) / 2
        for j in range(n1):
            temp = wn[j] - wn2[i]
            temp = -(temp * temp)
            res[i][j] = math.exp(temp / sigma)
    return res


def resample(srf, input, wn, FWH, col):
    srfsum = np.sum(srf, 1)
    inputpool = srfsum < 0.000001
    srfsum = srfsum[:, np.newaxis]
    weight = srf / srfsum
    weight[inputpool, :] = 0
    ninput = np.size(wn, 0)
    ntran = np.size(FWH, 0)
    transpec = np.zeros([ntran, col])

    for i in range(col):
        spec = input[:, i]
        spec = spec[:, np.newaxis]
        ref = weight.T * spec
        tran = np.sum(ref, 0)
        # tran = tran[:, np.newaxis]

        transpec[:, i] = tran

        # for j in range(ntran):
        #     sum1 = 0
        #     for k in range(ninput):
        #         a = wn[k]
        #         b = FWH[j, 0] + FWH[j, 1]
        #         c = FWH[j, 0] - FWH[j, 1]
        #         if (wn[k] <= FWH[j, 0] + FWH[j, 1]) and (wn[j] >= FWH[j, 0] - FWH[j, 1]):
        #             sum1 += ref[k][j]
        #     transpec[j][i] = sum1

    return transpec


def transform_Image_Data(inputfile, FWHfile, savefile, inputwarpfile, isspa=False):
    if isspa != False:
        spam = float(isspa)
        inputwarp = gdal.Warp(inputwarpfile, inputfile, format='ENVI', xRes=spam, yRes=spam,
                              resampleAlg=gdal.GRA_Bilinear, outputType=gdal.GDT_Float32)
    else:
        inputwarp = gdal.Open(inputfile)
    wn = np.zeros(inputwarp.RasterCount)
    im_geotrans = inputwarp.GetGeoTransform()  # 仿射矩阵
    im_proj = inputwarp.GetProjection()  # 地图投影信息
    arrayinput = inputwarp.ReadAsArray()
    for i in range(1, inputwarp.RasterCount + 1):
        wn[i - 1] = float(inputwarp.GetRasterBand(i).GetMetadataItem("wavelength"))

    FWH = np.loadtxt(FWHfile, comments=['Column', 'ENVI'])
    srf = calsrf(wn, FWH)
    bands = np.size(FWH, 0)
    col = np.size(arrayinput, 1)
    row = np.size(arrayinput, 2)
    output = np.zeros([bands, col, row])
    for i in range(row):
        input1 = arrayinput[:, :, i]
        output[:, :, i] = resample(srf, input1, wn, FWH, col)
    driver = gdal.GetDriverByName("ENVI")
    # filename = savefile + '.img'
    outdata = driver.Create((savefile + '.img'), row, col, bands, gdal.GDT_Float32)
    outdata.SetGeoTransform(im_geotrans)  # 写入仿射变换参数
    outdata.SetProjection(im_proj)  # 写入投影
    hdrtxt = 'wavelength units = Nanometers\n'
    hdrtxt = hdrtxt + 'wavelength = {\n'
    for i in range(bands):
        outdata.GetRasterBand(i + 1).WriteArray(output[i])

        if (i + 1) % 5 == 0:
            hdrtxt = hdrtxt + str(FWH[i, 0]) + ',\n'
        else:
            hdrtxt = hdrtxt + str(FWH[i, 0]) + ','
    del outdata
    if bands % 5 == 0:
        hdrtxt = hdrtxt[:-2] + '}'
    else:
        hdrtxt = hdrtxt[:-1] + '}'
    savefilehdr = savefile + '.hdr'
    with open(savefilehdr, 'a+') as sfh:
        sfh.write(hdrtxt)  # 在hdr文件中加入一行关于波长的信息


isspa = 30  # 输入空间分辨率 or False
inputfile = r'C:/Users/agrs/Desktop/to xx/test/数据转换/图像数据/hymap15m2.img'
FWHfile = r'C:/Users/agrs/Desktop/to xx/py文件/test/数据转换/图像数据/中心波长半高宽.txt'
savefile = r'C:/Users/agrs/Desktop/to xx/py文件/test/数据转换/图像数据/ccc'
# 无用的中间变量（缓存），待优化。。。。。。。
inputwarpfile = r'C:\inputwarp.img'

transform_Image_Data(inputfile, FWHfile, savefile, inputwarpfile, isspa=False)
