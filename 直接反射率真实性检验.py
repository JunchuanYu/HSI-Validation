from osgeo import gdal, osr
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


def spectruecal(inputfile, truefile, geofile, savefile, wstart, wend, iscoff=False, issam=False, isdis=False):
    # inputfile = 'D:\ZSXJY\\test\反射率产品真实性检验\间接\gf5_279band'
    # truefile = 'D:\ZSXJY\\test\反射率产品真实性检验\直接\hymap_to_gf5_gp22.txt'
    # geofile = 'D:\ZSXJY\\test\反射率产品真实性检验\直接\geo.txt'
    # savefile = 'D:\ZSXJY\\test\反射率产品真实性检验\直接\c'
    # 文件路径不能为空，为空则停止，且显示缺少

    name = []
    x_y = []
    with open(geofile, "r") as f:
        for line in f.readlines():
            line = line.strip('\n')  # 去掉列表中每一个元素的换行符
            if line[0] == '#':
                continue
            else:
                list1 = line.split(' ')
                name.append(list1[0])
                x_y.append([float(list1[1]), float(list1[2])])
    location = np.array(x_y)
    # 读取为np数组
    inputf = gdal.Open(inputfile)
    array_input = inputf.ReadAsArray()

    gtf = inputf.GetGeoTransform()
    srs = osr.SpatialReference()
    proj = inputf.GetProjection()  # 地图投影信息
    srs.ImportFromWkt(proj)
    srs_lon_lat = srs.CloneGeogCS()
    ct = osr.CoordinateTransformation(srs_lon_lat, srs)

    # x_range = range(0, inputf.RasterXSize)
    # y_range = range(0, inputf.RasterYSize)
    # x, y = np.meshgrid(x_range, y_range)
    # lon = gtf[0] + x * gtf[1] + y * gtf[2]
    # lat = gtf[3] + x * gtf[4] + y * gtf[5]

    # input = np.loadtxt(inputfile,  comments=['Column', 'ENVI'])
    true = np.loadtxt(truefile, comments=['Column', 'ENVI'])
    wn = true[:, 0]

    # 读取地理信息

    # 截取波长范围的数组
    wstart = float(wstart)
    wend = float(wend)
    if (wstart <= 0) or (wstart == ''):
        wstart = wn[0]
    if (wend <= wstart) or (wend == ''):
        wend = wn[-1]

    wn[wn > wend] = -1
    inputsbool = wn >= wstart

    input = array_input[inputsbool, :, :]
    true = true[inputsbool, :]
    row = np.size(input, 0)  # 计算 X 的行数
    # col = np.size(input, 2)  # 计算 X 的列数

    # 读取地理信息
    npoint = np.size(location, 0)
    input2 = np.zeros((row, npoint + 1))
    # wn1 = true[:,0]
    # wn1 = wn1[:, np.newaxis]
    input2[:, 0] = true[:, 0]
    for i in range(npoint):
        new = ct.TransformPoint(location[i, 1], location[i, 0])
        x = round((new[0] - gtf[0]) / gtf[1])
        y = round((new[1] - gtf[3]) / gtf[5])
        input2[:, i + 1] = input[:, y, x]
    col = np.size(input2, 1)

    # 三个复选框
    coff = np.zeros(col - 1)
    sam = np.zeros(col - 1)
    dis = np.zeros(col - 1)

    if iscoff:
        coff = speccoff(input2, true, row, col)
    if issam:
        sam = specsam(input2, true, row, col)
    if isdis:
        dis = specdis(input2, true, row, col)
    # coff = coff[:, np.newaxis]
    # coff = np.reshape(coff,(-1,5))
    # sam = np.reshape(sam, (-1, 5))
    # dis = np.reshape(dis, (-1, 5))

    savetxtfile = savefile + '.txt'
    with open(savetxtfile, 'a+') as f:
        np.savetxt(f, coff, fmt='%.6f', delimiter=' ', newline='\n', comments='', header='相关系数：', encoding='utf-8')
        np.savetxt(f, sam, fmt='%.6f', delimiter=' ', newline='\n', comments='', header='光谱角：', encoding='utf-8')
        np.savetxt(f, dis, fmt='%.6f', delimiter=' ', newline='\n', comments='', header='欧氏距离：', encoding='utf-8')


# 文件路径
inputfile = r'C:\Users\agrs\Desktop\to xx\py文件\test\反射率产品真实性检验\直接\gf5_279band'
truefile = r'C:\Users\agrs\Desktop\to xx\py文件\test\反射率产品真实性检验\直接\hymap_to_gf5_gp22.txt'
geofile = r'C:\Users\agrs\Desktop\to xx\py文件\test\反射率产品真实性检验\直接\geo.txt'
savefile = r'C:\Users\agrs\Desktop\to xx\py文件\test\反射率产品真实性检验\直接\ccc'

# 输入检验波长范围，wstart ~ wend（以nm为单位）
wstart = 800
wend = 1000

# 三个评价指标，默认值=Fales（光谱相关系数，光谱角，欧式距离）
iscoff = True
issam = True
isdis = True


spectruecal(inputfile, truefile, geofile, savefile, wstart, wend, iscoff, issam, isdis)
