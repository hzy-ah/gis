###### AVHRR
import numpy as np
import pandas as pd
from Vegetation_Index import *


def cloud_mask(reflectance1, reflectance2, reflectance3a):
    # 定义阈值
    threshold1 = 0.2
    threshold2 = 0.19
    threshold3 = 0.67
    threshold4 = 0.84
    threshold5 = 0.16
    threshold6 = 0.31
    threshold7 = 0.24
    threshold8 = 0.69

    # 创建空的云掩膜
    isCloud = np.zeros_like(reflectance1, dtype=np.uint8)

    # 根据阈值进行云检测
    isCloud[(reflectance1 > threshold1) & (reflectance2 > threshold2) & (reflectance3a < threshold3)] = 1
    isCloud[(reflectance1 > threshold1) & (reflectance2 > threshold5) & (reflectance3a > threshold4)] = 1
    isCloud[(reflectance1 > threshold6) & (reflectance2 > threshold7) & (reflectance3a > threshold8)] = 1

    return isCloud


def fmask(red, nir, swir1, swir2, bt_tir1, ndvi, ndsi):
    ###### AVHRR ######
    # red     1
    # nir     2
    # blue    1
    # green   1
    # swir1   3a
    # swir2   3b???
    # bt_tir1 bt4
    # 以下亮温单位为摄氏度，转换：单位K亮温 - 273.15
    bt_tir1 = bt_tir1 - 273.15
    # return cloud 1(ture)为云
    green = red
    blue = red

    olderr = np.seterr(all='ignore')

    ###### step 1 ######
    ### Basic Test
    basic_Test = np.logical_and(np.logical_and(swir2 > 0.03, bt_tir1 < 27),
                                np.logical_and(ndsi < 0.8, ndvi < 0.8))

    # ### Whiteness Test
    # meanvis = (red + blue + green) / 3
    # whiteness = np.zeros(red.shape)
    # for i in [red, blue, green]:
    #     whiteness0 = abs((i - meanvis) / meanvis)
    #     whiteness = whiteness + whiteness0
    # whiteness_Test = whiteness < 0.7

    # ### HOT Test
    # hot_Test = (blue - 0.5 * red - 0.08) > 0

    ### Ratio Test
    ratio_Test = (nir / swir1) > 0.75

    ### step 1 pass pcp
    # pcp = np.logical_and(np.logical_and(basic_Test, whiteness_Test),
    #                      np.logical_and(hot_Test, ratio_Test))
    pcp = np.logical_and(basic_Test, ratio_Test)

    ###### step 2 ######
    ### WATER
    water_Test = np.logical_or(np.logical_and(ndvi < 0.01, nir < 0.11),
                               np.logical_and(ndvi < 0.1, nir < 0.05))
    water_Test_1 = water_Test < 1
    clear_sky_water = np.logical_and(water_Test, swir2 < 0.03)
    bt_tir1_water = bt_tir1[clear_sky_water]  # clear
    bt_tir1_water1 = pd.Series(bt_tir1_water)
    T_water = bt_tir1_water1.quantile(0.825)
    wTemp_prob = (T_water - bt_tir1) / 4
    wBri_prob = np.minimum(swir1, 0.11) / 0.11
    wCloud_prob = wTemp_prob * wBri_prob

    pcp1 = wCloud_prob > 0.5  ### 0.5 可调试

    ### LAND
    clear_sky_land = np.logical_and(pcp < 1, water_Test_1)
    bt_tir1_land = bt_tir1[clear_sky_land]
    bt_tir1_land1 = pd.Series(bt_tir1_land)
    T_high = bt_tir1_land1.quantile(0.825)
    T_low = bt_tir1_land1.quantile(0.175)
    # T_high = np.percentile(bt_tir1_land, 82.5)
    # T_low = np.percentile(bt_tir1_land, 17.5)
    lTemp_prob = (T_high + 4 - bt_tir1) / (T_high + 4 - (T_low - 4))
    # lVar_prob = 1 - np.maximum(abs(ndvi), abs(ndsi), whiteness)
    lVar_prob = 1 - np.maximum(abs(ndvi), abs(ndsi))
    lCloud_prob = lTemp_prob * lVar_prob
    lCloud_prob_land = lCloud_prob[clear_sky_land]  # clear
    lCloud_prob_land1 = pd.Series(lCloud_prob_land)
    land_threshold = 0.2 + lCloud_prob_land1.quantile(0.825)  ### 0.2 可调试

    pcp2 = lCloud_prob > land_threshold

    cloud = np.logical_or(np.logical_or((pcp & water_Test & pcp1),
                                        (pcp & water_Test_1 & pcp2)),
                          np.logical_or(((lCloud_prob > 0.99) & water_Test_1),
                                        (bt_tir1 < (T_low - 35))))  ### 35 可调试

    return cloud


# 计算云掩膜
def Get_isCloud(TOA, BT):
    red = TOA[0]
    nir = TOA[1]
    blue = TOA[0]
    green = TOA[0]
    swir1 = TOA[2]
    swir2 = TOA[2]
    bt_tir1 = BT[0]
    ndvi = NDVI(nir, red)
    ndsi = NDSI(swir1, green)

    isCloud = cloud_mask(TOA[0], TOA[1], TOA[2])
    # isCloud = np.zeros((nir.shape[0], nir.shape[1]))
    # for i in range(red.shape[0]):
    #     for j in range(red.shape[1]):
    #         isCloud[i][j] = fmask(red[i][j], nir[i][j], swir1[i][j], swir2[i][j],
    #                               bt_tir1[i][j], ndvi[i][j], ndsi[i][j])

    return isCloud
