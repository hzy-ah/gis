import pandas as pd
import numpy as np
from osgeo import gdal
from Vegetation_Index import *


# MODIS
def CLOUD(red, nir, swir):
    # 设定阈值进行云检测
    cloud_mask = np.logical_and(nir > 0.35, np.logical_and(red > 0.11, swir > 0.08))
    return cloud_mask.astype(int)


# GK2A
def CLOUD_GK2A(bands):
    # 设置阈值
    thresholds = [0.32, 0.32, 0.32, 0.35, 0.3, 0.3]
    # 执行云检测
    cloud_mask = np.zeros_like(bands[0])  # 创建一个和输入文件大小相同的零矩阵

    for band, threshold in zip(bands, thresholds):
        # 使用阈值法进行云检测
        cloud_pixels = np.where(band > threshold, 1, 0)
        cloud_mask += cloud_pixels
    cloud_mask = np.where(cloud_mask > 0, 1, 0)
    return cloud_mask


def fmask(red, nir, blue, green, swir1, swir2, cirrus, bt_tir1, ndvi, ndsi):
    ### AQUA MODIS 1.64|可能有问题，需采用插值或其他方法获取，目前照旧使用
    ######  MODIS VIIRS FY3D_MERSI FY3C_VIRR GK2A_AMI ######
    # red     1     I1      3          1         3
    # nir     2     I2      4          2         4
    # blue    3     M2      1          7         1
    # green   4     M4      2          9         2
    # swir1   6     I3      6          6         6
    # swir2   7     M11     7          3???      7
    # cirrus  26    M9      5          10        5
    # bt_tir1 bt31  btM15   bt24       bt4       bt14
    # 以下亮温单位为摄氏度，转换：单位K亮温 - 273.15
    bt_tir1 = bt_tir1 - 273.15
    # return cloud 1(ture)为云

    olderr = np.seterr(all='ignore')

    ###### step 1 ######
    ### Basic Test
    cirrus_prob = cirrus / 0.04
    basic_Test = np.logical_and(np.logical_and(swir2 > 0.03, bt_tir1 < 27),
                                np.logical_and(ndsi < 0.8, ndvi < 0.8))

    ### Whiteness Test
    meanvis = (red + blue + green) / 3
    whiteness = np.zeros(red.shape)
    for i in [red, blue, green]:
        whiteness0 = abs((i - meanvis) / meanvis)
        whiteness = whiteness + whiteness0
    whiteness_Test = whiteness < 0.7

    ### HOT Test
    hot_Test = (blue - 0.5 * red - 0.08) > 0

    ### Ratio Test
    ratio_Test = (nir / swir1) > 0.75

    ### step 1 pass pcp
    pcp = np.logical_and(np.logical_and(basic_Test, whiteness_Test),
                         np.logical_and(hot_Test, ratio_Test))

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
    wCloud_prob = wTemp_prob * wBri_prob + cirrus_prob

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
    lVar_prob = 1 - max(abs(ndvi), abs(ndsi), whiteness)

    lCloud_prob = lTemp_prob * lVar_prob + cirrus_prob
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
def Get_isCloud(TOA_500, TOA_500_26, BT_500_31):
    red = TOA_500[0]
    nir = TOA_500[1]
    blue = TOA_500[2]
    green = TOA_500[3]
    swir1 = TOA_500[5]
    swir2 = TOA_500[6]
    cirrus = TOA_500_26
    bt_tir1 = BT_500_31
    ndvi = NDVI(nir, red)
    ndsi = NDSI(swir1, green)

    isCloud = CLOUD(red, nir, swir1)

    return isCloud


def Get_isCloud_AMI(TOA1, TOA2):
    bands = []
    for i in range(len(TOA1)):
        bands.append(TOA1[i])
    for j in range(len(TOA2)):
        bands.append(TOA2[j])

    isCloud = CLOUD_GK2A(bands)
    return isCloud

