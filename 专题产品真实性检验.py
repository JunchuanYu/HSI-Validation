import numpy as np
import math
from sklearn.metrics import confusion_matrix
from osgeo import gdal


def compute_kappa(array_input, array_true):
    nclass = np.unique(array_true)  # 专题文件中有哪些类别
    n = np.size(nclass)  # 专题文件中有几种类别

    # 为防止裁剪过程两个数据大小不一致
    sample1 = np.size(array_true, 1)
    line1 = np.size(array_true, 0)
    sample2 = np.size(array_input, 1)
    line2 = np.size(array_input, 0)
    if (sample1 == sample2) and (line1 == line2):
        sample = sample1
        line = line1
    else:
        sample = min(sample1, sample2)
        line = min(line1, line2)
        array_true = array_true[0:line, 0:sample]
        array_input = array_input[0:line, 0:sample]

    # 将二维的数组变成一列
    array_true1 = array_true.flatten()
    array_input1 = array_input.flatten()
    # 混淆矩阵计算
    cm = confusion_matrix(array_true1, array_input1)  # cm为混淆矩阵

    sumlie = np.sum(cm, 0)  # 混淆矩阵列和
    sumhang = np.sum(cm, 1)  # 混淆矩阵行和
    tri = cm[np.arange(n), np.arange(n)]  # 矩阵元素的对角线提取

    OA = np.sum(tri) / (sample * line)  # 总体精度
    pe = np.sum(sumlie * sumhang) / ((sample * line) ** 2)
    kappa = (OA - pe) / (1 - pe)  # kappa系数
    prediction = np.zeros(n)  # 预测精度
    recall = np.zeros(n)  # 召回率
    f1 = np.zeros(n)
    for i in range(n):
        i_pre_sum = cm.sum(axis=0)[i]
        i_correct_sum = cm[i][i]
        i_recall_sum = cm.sum(axis=1)[i]
        prediction[i] = 0
        if i_pre_sum != 0:
            prediction[i] = 100 * float(i_correct_sum) / float(i_pre_sum)
        recall[i] = 0
        if i_recall_sum != 0:
            recall[i] = 100 * float(i_correct_sum) / float(i_recall_sum)
        if (prediction[i] + recall[i]) == 0:
            f1[i] = 0
        f1[i] = 2 * prediction[i] * recall[i] / (prediction[i] + recall[i])
    # p = prediction.sum() / len(prediction)
    # r = recall.sum() / len(recall)

    return n, nclass, cm, OA, kappa, prediction, recall, f1


def classtrue(inputfile, truefile, savefile, shpfile=False, inputwarpfile1=False, inputwarpfile2=False):
    gdal.AllRegister()

    # 文件路径不能为空，为空则停止，且显示缺少
    if not shpfile:
        tshp = 0
    else:
        tshp = 1
        # if tshp == 1:
        #     trueshpfile = shpfile

    input = gdal.Open(inputfile)
    true = gdal.Open(truefile)
    # shp = shapefile.Reader(shpfile)
    # trueshp = shapefile.Reader(trueshpfile)
    if tshp == 1:
        inputclass = gdal.Warp(inputwarpfile1, input, format='GTiff', cutlineDSName=shpfile, cropToCutline=True,
                               dstNodata=0)
        trueclass = gdal.Warp(inputwarpfile2, true, format='GTiff', cutlineDSName=shpfile, cropToCutline=True,
                              dstNodata=0)
        array_input = inputclass.GetRasterBand(1).ReadAsArray()
        array_true = trueclass.GetRasterBand(1).ReadAsArray()
    else:
        array_input = input.GetRasterBand(1).ReadAsArray()
        array_true = true.GetRasterBand(1).ReadAsArray()

    sample1 = np.size(array_true, 1)
    line1 = np.size(array_true, 0)
    sample2 = np.size(array_input, 1)
    line2 = np.size(array_input, 0)
    if (sample1 == sample2) and (line1 == line2):
        sample = sample1
        line = line1
    else:
        sample = min(sample1, sample2)
        line = min(line1, line2)
        array_true = array_true[0:line, 0:sample]
        array_input = array_input[0:line, 0:sample]

    # 建立混淆矩阵
    n, nclass, cm, OA, kappa, prediction, recall, f1 = compute_kappa(array_input, array_true)

    OA = round(OA, 6)
    kappa = round(kappa, 6)

    savetxtfile = savefile + '.txt'
    with open(savetxtfile, 'a+') as f:
        writetxt = '总体分类精度为：' + str(OA) + '\n' + 'kappa系数为：' + str(kappa) + '\n'
        f.write(writetxt)
        np.savetxt(f, cm, fmt='%d', delimiter='\t', newline='\n', comments='', header='混淆矩阵：', encoding='utf-8')
        writetxt1 = '类别\t' + '对应图中值\t' + '精确度\t' + '召回率\t' + 'f1\n'
        f.write(writetxt1)
        for i in range(n):
            prediction[i] = round(prediction[i], 2)
            recall[i] = round(recall[i], 2)
            f1[i] = round(f1[i], 2)
            writetxtclass = str(i + 1) + '\t' + str(nclass[i]) + '\t' + str(prediction[i]) + '\t' + str(
                recall[i]) + '\t' + str(f1[i]) + '\n'
            f.write(writetxtclass)
        # np.savetxt(f, sam, fmt='%.6f', delimiter=' ', newline='\n', comments='', header='光谱角：', encoding='utf-8')
        # np.savetxt(f, dis, fmt='%.6f', delimiter=' ', newline='\n', comments='', header='欧氏距离：', encoding='utf-8')
    print("OA: ", OA)
    print("kappa: ", kappa)


# 文件路径
inputfile = r'C:/Users/agrs/Desktop/to xx/test/专题产品真实性检验/10mcankaolvhuanglan.img'
truefile = r'C:/Users/agrs/Desktop/to xx/test/专题产品真实性检验/15m_10nm_分类图_lvhuanglan.img'
# 中文目录可能会报错
savefile = r'C:/Users/agrs/Desktop/to xx/test/专题产品真实性检验/ccc'
# shpfile 默认值=False，输入待真实性检验专题产品SHP文件，传入地址，对应 app 中的 “检测范围”
shpfile = r'C:/Users/agrs/Desktop/to xx/test/专题产品真实性检验/ZTCP.shp'
# inputwarp, 缓存文件
inputwarpfile1 = r'C:/Users/agrs/Desktop/to xx/example/warp/inputwarp1.tif'
inputwarpfile2 = r'C:/Users/agrs/Desktop/to xx/example/warp/inputwarp2.tif'


classtrue(inputfile, truefile, savefile, shpfile, inputwarpfile1, inputwarpfile2)
