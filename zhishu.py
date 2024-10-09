import numpy as np
import math


def RVI(NIR, R):
    rvi = np.zeros((NIR.shape[0], NIR.shape[1]))
    for i in range(NIR.shape[0]):
        for j in range(NIR.shape[1]):
            if R[i][j] == -999 or R[i][j] == 0 or NIR[i][j] == -999:
                rvi[i][j] = -999
            else:
                rvi[i][j] = NIR[i][j] / R[i][j]

    return rvi


def NDVI(NIR, R):
    ndvi = np.zeros((NIR.shape[0], NIR.shape[1]))
    for i in range(NIR.shape[0]):
        for j in range(NIR.shape[1]):
            if NIR[i][j] + R[i][j] == 0 or R[i][j] == -999 or NIR[i][j] == -999:
                ndvi[i][j] = -999
            else:
                ndvi[i][j] = -(R[i][j] - NIR[i][j]) / (NIR[i][j] + R[i][j])
                if ndvi[i][j] > 1:
                    ndvi[i][j] = 1

    return ndvi


def NDSI(SWIR, GREEN):
    ndsi = np.zeros((SWIR.shape[0], SWIR.shape[1]))
    for i in range(SWIR.shape[0]):
        for j in range(SWIR.shape[1]):
            if (SWIR[i][j] + GREEN[i][j]) == 0 or SWIR[i][j] == -999 or GREEN[i][j] == -999:
                ndsi[i][j] = -999
            else:
                ndsi[i][j] = (GREEN[i][j] - SWIR[i][j]) / (SWIR[i][j] + GREEN[i][j])


    return ndsi


def SAVI(NIR, R, L):
    savi = - 999 * np.ones((NIR.shape[0], NIR.shape[1]))
    for i in range(NIR.shape[0]):
        for j in range(NIR.shape[1]):
            if R[i][j] == -999 or NIR[i][j] < 0 or R[i][j] < 0 or NIR[i][j] > 1 or R[i][j] > 1:
                savi[i][j] = -999
            else:
                savi[i][j] = (NIR[i][j] - R[i][j]) * (1 + L) / (NIR[i][j] + R[i][j] + L)



    return savi


def GVI(Band1, Band2, Band3, Band4):
    gvi = np.zeros((Band1.shape[0], Band1.shape[1]))
    for i in range(Band1.shape[0]):
        for j in range(Band1.shape[1]):
            if Band1[i][j] < 0 or Band2[i][j] < 0 or Band3[i][j] < 0 or Band4[i][j] < 0 or \
                    Band1[i][j] > 1 or Band2[i][j] > 1 or Band3[i][j] > 1 or Band4[i][j] > 1:
                gvi[i][j] = -999
            else:
                gvi[i][j] = -0.283 * Band1[i][j] - 0.66 * Band2[i][j] + 0.577 * Band3[i][j] + 0.388 * Band4[i][j]

    return gvi


def EVI(NIR, RED, BLUE):
    evi = np.zeros((NIR.shape[0], NIR.shape[1]))
    for i in range(NIR.shape[0]):
        for j in range(NIR.shape[1]):
            if NIR[i][j] < 0 or RED[i][j] < 0 or BLUE[i][j] < 0 or NIR[i][j] > 1 or RED[i][j] > 1 or BLUE[i][j] > 1:
                evi[i][j] = -999
            else:
                evi[i][j] = 2.5 * (NIR[i][j] - RED[i][j]) / (NIR[i][j] + 6 * RED[i][j] - 7.5 * BLUE[i][j] - 1)


    return evi


def LAI(NDVI, Landcover):
    result = np.zeros((NDVI.shape[0], NDVI.shape[1]))
    for i in range(NDVI.shape[0]):
        for j in range(NDVI.shape[1]):

            if Landcover[i][j] == 0:
                result[i][j] = -999

            if Landcover[i][j] == 1:
                if NDVI[i][j] < 0.125:
                    result[i][j] = -999
                if 0.125 <= NDVI[i][j] <= 0.825:
                    result[i][j] = 0.1836 * math.exp(4.37 * NDVI[i][j])
                if 0.825 < NDVI[i][j]:
                    result[i][j] = 6.606

            if Landcover[i][j] == 2:
                if NDVI[i][j] < 0.125:
                    result[i][j] = -999
                if 0.125 <= NDVI[i][j] <= 0.825:
                    result[i][j] = 0.0884 * math.exp(4.96 * NDVI[i][j])
                if 0.825 < NDVI[i][j]:
                    result[i][j] = 6.091

    return result
