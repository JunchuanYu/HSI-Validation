import numpy as np
import math


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


#  光谱数据转换
def transform_Spectral_Data(inputfile, FWHfile, savefile):
    FWH = np.loadtxt(FWHfile, comments=['Column', 'ENVI'])
    input = np.loadtxt(inputfile, comments=['Column', 'ENVI'])
    wn = input[:, 0]
    input = input[:, 1:]
    srf = calsrf(wn, FWH)
    col = np.size(input, 1)
    output = resample(srf, input, wn, FWH, col)
    savetxtfile = savefile + '.txt'
    # open(savetxtfile, 'ab') as f
    wn2 = FWH[:, 0]
    wn2 = wn2[:, np.newaxis]
    output1 = np.hstack((wn2, output))
    np.savetxt(savetxtfile, output1, fmt='%.6f', delimiter='\t', newline='\n', comments='', encoding='utf-8')


inputfile = r'C:\Users\agrs\Desktop\to xx\py文件\test\数据转换\光谱数据\hymap_gp2.txt'
FWHfile = r'C:\Users\agrs\Desktop\to xx\py文件\test\数据转换\光谱数据\中心波长半高宽.txt'
savefile = r'C:\Users\agrs\Desktop\to xx\py文件\test\数据转换\光谱数据\\ccc'

transform_Spectral_Data(inputfile, FWHfile, savefile)
