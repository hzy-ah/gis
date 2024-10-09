import numpy as np
import netCDF4 as nc
import math


# 比值植被指数
# NIR为近红外波段反射率，R是红光波段反射率
def RVI(NIR, R):
    rvi = -999 * np.ones((NIR.shape[0], NIR.shape[1]))
    for i in range(NIR.shape[0]):
        for j in range(NIR.shape[1]):
            if R[i][j] != 0 and R[i][j] != -999 and NIR[i][j] != -999:
                rvi[i][j] = NIR[i][j] / R[i][j]
    return rvi


# 归一化植被指数
def NDVI(NIR, R):
    ndvi = -999 * np.ones((NIR.shape[0], NIR.shape[1]))
    for i in range(NIR.shape[0]):
        for j in range(NIR.shape[1]):
            if NIR[i][j] + R[i][j] != 0 and R[i][j] != -999 and NIR[i][j] != -999:
                ndvi[i][j] = (NIR[i][j] - R[i][j]) / (NIR[i][j] + R[i][j])

    return ndvi


# 土壤调整植被指数
# L=0 时，表示植被覆盖度为零；L=1时，表示土壤背景的影响为零，即植被覆盖度非常高，土壤背景的影响为零，这种情况只有在被树冠浓密的高大树木覆盖的地方才会出现。
def SAVI(NIR, R, L):
    savi = -999 * np.ones((NIR.shape[0], NIR.shape[1]))
    for i in range(NIR.shape[0]):
        for j in range(NIR.shape[1]):
            if NIR[i][j] + R[i][j] + L != 0 and R[i][j] != -999 and NIR[i][j] != -999:
                savi[i][j] = (NIR[i][j] - R[i][j]) * (1 + L) / (NIR[i][j] + R[i][j] + L)

    return savi


# 绿度植被指数-0.2848 * Band1 - 0.2435 * Band2 - 0.5436 * Band3 + 0.7243 * Band4 + 0.084 * Band5 - 1.18 * Band7
# band1蓝绿波段，band2红波段，band3，band4红外波段
def GVI(Band1, Band2, Band3, Band4):
    gvi = -999 * np.ones((Band1.shape[0], Band1.shape[1]))
    for i in range(Band1.shape[0]):
        for j in range(Band1.shape[1]):
            if Band1[i][j] != -999 and Band2[i][j] != -999 and Band3[i][j] != -999 and Band4[i][j] != -999:
                gvi[i][j] = -0.283 * Band1[i][j] - 0.66 * Band2[i][j] + 0.577 * Band3[i][j] + 0.388 * Band4[i][j]
    return gvi


# 增强植被指数
# NIR,RED,BLUE分别代表近红外波段、红光波段、蓝光波段的反射率
def EVI(NIR, RED, BLUE):
    evi = -999 * np.ones((NIR.shape[0], NIR.shape[1]))
    for i in range(NIR.shape[0]):
        for j in range(NIR.shape[1]):
            if NIR[i][j] + 6 * RED[i][j] - 7.5 * BLUE[i][j] - 1 != 0 and NIR[i][j] != -999 and RED[i][j] != -999 and \
                    BLUE[i][j] != -999:
                evi[i][j] = 2.5 * (NIR[i][j] - RED[i][j]) / (NIR[i][j] + 6 * RED[i][j] - 7.5 * BLUE[i][j] - 1)

    return evi


# 基于经验算法估算叶面积指数
def LAI(NDVI, Landcover):
    result = -999 * np.ones((NDVI.shape[0], NDVI.shape[1]))
    for i in range(NDVI.shape[0]):
        for j in range(NDVI.shape[1]):

            if NDVI[i][j] != -999:
                if Landcover[i][j] == 0:  # 水体
                    result[i][j] = 0

                if Landcover[i][j] == 1:  # 林地
                    if NDVI[i][j] < 0.125:
                        result[i][j] = 0
                    if 0.125 <= NDVI[i][j] <= 0.825:
                        result[i][j] = 0.1836 * math.exp(4.37 * NDVI[i][j])
                    if 0.825 < NDVI[i][j]:
                        result[i][j] = 6.606

                if Landcover[i][j] == 2:  # 草地
                    if NDVI[i][j] < 0.125:
                        result[i][j] = 0
                    if 0.125 <= NDVI[i][j] <= 0.825:
                        result[i][j] = 0.0884 * math.exp(4.96 * NDVI[i][j])
                    if 0.825 < NDVI[i][j]:
                        result[i][j] = 6.091
    return result


# 归一化差雪指数（云检测从算法需要）
def NDSI(SWIR, GREEN):
    ndsi = -999 * np.ones((SWIR.shape[0], SWIR.shape[1]))
    for i in range(SWIR.shape[0]):
        for j in range(SWIR.shape[1]):
            if SWIR[i][j] + GREEN[i][j] != 0 and SWIR[i][j] != -999 and GREEN[i][j] != -999:
                ndsi[i][j] = (GREEN[i][j] - SWIR[i][j]) / (SWIR[i][j] + GREEN[i][j])

    return ndsi


def RVI_nc(rvi, lat, lon, Ins_Dir, date):
    Output_path = Ins_Dir + '/' + 'RVI' + date + '.nc'
    nc_file = nc.Dataset(Output_path, 'w', format='NETCDF4')
    nc_file.createDimension('lat', len(lat))
    nc_file.createDimension('lon', len(lon))
    nc_file.createDimension('size', 1)

    group = nc_file.createGroup('RVI')

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon

    RVI = group.createVariable('RVI', np.float32,
                               ('size', 'lat', 'lon'))
    group.variables['RVI'][:] = rvi
    RVI.time = date
    RVI.band_names = ['RVI']


def NDVI_nc(ndvi, lat, lon, Ins_Dir, date):
    Output_path = Ins_Dir + '/' + 'NDVI' + date + '.nc'
    nc_file = nc.Dataset(Output_path, 'w', format='NETCDF4')
    nc_file.createDimension('lat', len(lat))
    nc_file.createDimension('lon', len(lon))
    nc_file.createDimension('size', 1)

    group = nc_file.createGroup('NDVI')

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon

    NDVI = group.createVariable('NDVI', np.float32,
                                ('size', 'lat', 'lon'))
    group.variables['NDVI'][:] = ndvi
    NDVI.time = date
    NDVI.band_names = ['NDVI']


def SAVI_nc(savi, lat, lon, Ins_Dir, date):
    Output_path = Ins_Dir + '/' + 'SAVI' + date + '.nc'
    nc_file = nc.Dataset(Output_path, 'w', format='NETCDF4')
    nc_file.createDimension('lat', len(lat))
    nc_file.createDimension('lon', len(lon))
    nc_file.createDimension('size', 1)

    group = nc_file.createGroup('SAVI')

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon

    SAVI = group.createVariable('SAVI', np.float32,
                                ('size', 'lat', 'lon'))
    group.variables['SAVI'][:] = savi
    SAVI.time = date
    SAVI.band_names = ['SAVI']


def GVI_nc(gvi, lat, lon, Ins_Dir, date):
    Output_path = Ins_Dir + '/' + 'GVI' + date + '.nc'
    nc_file = nc.Dataset(Output_path, 'w', format='NETCDF4')
    nc_file.createDimension('lat', len(lat))
    nc_file.createDimension('lon', len(lon))
    nc_file.createDimension('size', 1)

    group = nc_file.createGroup('GVI')

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon

    GVI = group.createVariable('GVI', np.float32,
                               ('size', 'lat', 'lon'))
    group.variables['GVI'][:] = gvi
    GVI.time = date
    GVI.band_names = ['GVI']


def EVI_nc(evi, lat, lon, Ins_Dir, date):
    Output_path = Ins_Dir + '/' + 'EVI' + date + '.nc'
    nc_file = nc.Dataset(Output_path, 'w', format='NETCDF4')
    nc_file.createDimension('lat', len(lat))
    nc_file.createDimension('lon', len(lon))
    nc_file.createDimension('size', 1)

    group = nc_file.createGroup('EVI')

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon

    EVI = group.createVariable('EVI', np.float32,
                               ('size', 'lat', 'lon'))
    group.variables['EVI'][:] = evi
    EVI.time = date
    EVI.band_names = ['EVI']


def LAI_nc(lai, lat, lon, Ins_Dir, date):
    Output_path = Ins_Dir + '/' + 'LAI' + date + '.nc'
    nc_file = nc.Dataset(Output_path, 'w', format='NETCDF4')
    nc_file.createDimension('lat', len(lat))
    nc_file.createDimension('lon', len(lon))
    nc_file.createDimension('size', 1)

    group = nc_file.createGroup('LAI')

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon

    LAI = group.createVariable('LAI', np.float32,
                               ('size', 'lat', 'lon'))
    group.variables['LAI'][:] = lai
    LAI.time = date
    LAI.band_names = ['LAI']
