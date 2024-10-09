import numpy as np
import pandas as pd


def fmask(red, nir, blue, green, swir1, swir2, cirrus, bt_tir1, ndvi, ndsi):
    bt_tir1 = bt_tir1 - 273.15

    olderr = np.seterr(all='ignore')

    cirrus_prob = cirrus / 0.04
    basic_Test = np.logical_and(np.logical_and(swir2 > 0.03, bt_tir1 < 27),
                                np.logical_and(ndsi < 0.8, ndvi < 0.8))

    meanvis = (red + blue + green) / 3
    whiteness = np.zeros(red.shape)
    for i in [red, blue, green]:
        whiteness0 = abs((i - meanvis) / meanvis)
        whiteness = whiteness + whiteness0
    whiteness_Test = whiteness < 0.7

    hot_Test = (blue - 0.5 * red - 0.08) > 0

    ratio_Test = (nir / swir1) > 0.75

    pcp = np.logical_and(np.logical_and(basic_Test, whiteness_Test),
                         np.logical_and(hot_Test, ratio_Test))

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

    pcp1 = wCloud_prob > 0.5

    clear_sky_land = np.logical_and(pcp < 1, water_Test_1)
    bt_tir1_land = bt_tir1[clear_sky_land]
    bt_tir1_land1 = pd.Series(bt_tir1_land)
    T_high = bt_tir1_land1.quantile(0.825)
    T_low = bt_tir1_land1.quantile(0.175)
    lTemp_prob = (T_high + 4 - bt_tir1) / (T_high + 4 - (T_low - 4))
    lVar_prob = 1 - np.maximum(abs(ndvi), abs(ndsi), whiteness)
    lCloud_prob = lTemp_prob * lVar_prob + cirrus_prob
    lCloud_prob_land = lCloud_prob[clear_sky_land]
    lCloud_prob_land1 = pd.Series(lCloud_prob_land)
    land_threshold = 0.2 + lCloud_prob_land1.quantile(0.825)

    pcp2 = lCloud_prob > land_threshold

    cloud = np.logical_or(np.logical_or((pcp & water_Test & pcp1),
                                        (pcp & water_Test_1 & pcp2)),
                          np.logical_or(((lCloud_prob > 0.99) & water_Test_1),
                                        (bt_tir1 < (T_low - 35))))

    return cloud


def fmask2(nir_band, red_band, swir_band, ndvi):
    cloud_mask = np.logical_and(nir_band > 0.35, np.logical_and(red_band > 0.11, swir_band > 0.08))
    cloud_mask = np.logical_and(cloud_mask, np.logical_not(red_band == -999))
    cloud_mask = cloud_mask.astype(int)

    return cloud_mask
