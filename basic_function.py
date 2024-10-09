import numpy as np
from osgeo import gdal
import h5py
import os
import netCDF4 as nc
import read_data as rd
import zhishu as zs


# 判断数据是否有足够的像素点落在制定区域
def ifuse(geo_path, band_path, data_type):
    if data_type == 'FY3C':
        band = h5py.File(band_path, 'r')
        drn = band.attrs['Day Or Night Flag']
        if drn[0] != 68:
            return -1

        geo = h5py.File(geo_path, 'r')
        h5_geo_Lat = geo['Geolocation']['Latitude']
        h5_geo_Lon = geo['Geolocation']['Longitude']
        lat = np.array(h5_geo_Lat).astype(float)  # 纬度
        lon = np.array(h5_geo_Lon).astype(float)  # 经度
        get_piex = 0
        for i in range(lat.shape[0]):
            for j in range(lat.shape[1]):
                if 22 <= lat[i][j] <= 29 and 115 <= lon[i][j] <= 122:
                    get_piex = get_piex + 1
        if get_piex < 50:
            return -1

        return 0

    if data_type == 'FY3D':
        band = h5py.File(band_path, 'r')
        drn = band.attrs['Day Or Night Flag']
        if drn[0] != 68:
            return -1
        geo = h5py.File(geo_path, 'r')
        h5_geo_Lat = geo['Geolocation']['Latitude']
        h5_geo_Lon = geo['Geolocation']['Longitude']
        lat = np.array(h5_geo_Lat).astype(float)  # 纬度
        lon = np.array(h5_geo_Lon).astype(float)  # 经度

        get_piex = 0
        for i in range(lat.shape[0]):
            for j in range(lat.shape[1]):
                if 22 <= lat[i][j] <= 29 and 115 <= lon[i][j] <= 122:
                    get_piex = get_piex + 1
        if get_piex < 50:
            return -1

        return 0

    if data_type == 'VIIRS':
        geo = h5py.File(geo_path, 'r')
        allkey = []
        for key in geo.keys():
            allkey.append(key)

        allkey2 = []
        for key in geo[allkey[0]].keys():
            allkey2.append(key)

        h5_geo_Lat = geo[allkey[0]][allkey2[0]]['Latitude']
        h5_geo_Lon = geo[allkey[0]][allkey2[0]]['Longitude']

        lat = np.array(h5_geo_Lat).astype(float)  # 纬度
        lon = np.array(h5_geo_Lon).astype(float)  # 经度

        get_piex = 0
        for i in range(lat.shape[0]):
            for j in range(lat.shape[1]):
                if 22 <= lat[i][j] <= 29 and 115 <= lon[i][j] <= 122:
                    get_piex = get_piex + 1
        if get_piex < 50:
            return -1

        return 0


# 读取经纬度等数据
def geo(geo_filename, data_type):
    geo = h5py.File(geo_filename, 'r')
    if data_type == 'FY3D_250':
        lat, lon, data_range = rd.read_geo(geo, data_type)
        return lat, lon, data_range

    h, lat, lon, sa_a, sa_z, so_a, so_z = rd.read_geo(geo, data_type)

    A = rd.relative_angle(sa_a, sa_z, so_a, so_z)

    return h, lat, lon, sa_z, so_z, A


# 读取卫星波段数据
def band(band_path, data_type, band_type=0, data_range=None):
    data = h5py.File(band_path, 'r')
    band = rd.read_band(data, data_type, band_type, data_range=data_range)
    return band


# 读取查找表的有用信息
def table_data(filename):
    data = nc.Dataset(filename)

    Pa = data.variables['Pa']
    Ra = data.variables['Ra']
    TT = data.variables['TT']

    AOD = data.variables['AOD550nm']
    Altitude = data.variables['Altitude']
    Mu = data.variables['VZA']
    Mu0 = data.variables['SZA']
    RAA = data.variables['RAA']

    Pa = np.array(Pa)
    Ra = np.array(Ra)
    TT = np.array(TT)

    AOD = np.array(AOD)
    Altitude = np.array(Altitude)
    Mu = np.array(Mu)
    Mu0 = np.array(Mu0)
    RAA = np.array(RAA)

    return Pa, Ra, TT, AOD, Altitude, Mu, Mu0, RAA


# 几何校正函数
def swath_georeference_vrt(out_dir, name, outputbounds, xres, yres):
    out_filename = out_dir + '/' + 'after' + str(name) + '.tif'
    vrt_path = out_dir + '/' + 'vrt.vrt'
    ds = gdal.Warp(out_filename, vrt_path, format='GTiff', outputBounds=outputbounds,
                   dstSRS='EPSG:4326', xRes=xres, yRes=yres, geoloc=True,
                   dstNodata=-999, srcNodata=-32768,
                   resampleAlg=gdal.gdalconst.GRA_Bilinear)
    os.remove(vrt_path)

    return ds, out_filename


# 大气校正数据准备函数
def atmos(ref, sa_z, so_z, A, table_file, cloud):
    path = './ndata/dem.tif'
    f = gdal.Open(path)
    h = f.ReadAsArray()
    h = np.array(h).astype(float)

    if cloud.shape[0] == 700:
        h = h[::2, ::2]

    Pa, Ra, TT, AOD, Altitude, Mu, Mu0, RAA = table_data(table_file)

    Rs = atmos_correct(ref, h, sa_z, so_z, A, Pa, Ra, TT, cloud)

    return Rs


# 大气校正函数，根据各波段的查找表进行校正
def atmos_correct(data, my_Altitude, sa_z, so_z, A, Pa, Ra, TT, cloud):
    Pa_PList = np.zeros((sa_z.shape[0], sa_z.shape[1]))
    Ra_PList = np.zeros((sa_z.shape[0], sa_z.shape[1]))
    TT_PList = np.zeros((sa_z.shape[0], sa_z.shape[1]))
    Rs = np.zeros((sa_z.shape[0], sa_z.shape[1]))

    for i in range(sa_z.shape[0]):
        for j in range(sa_z.shape[1]):
            if cloud[i][j] == 0:
                if data[i][j] != -999:
                    H = my_Altitude[i][j]
                    my_Mu0 = int(so_z[i][j])
                    my_Mu = int(sa_z[i][j])
                    my_RAA = int(A[i][j])

                    if H > 3:
                        H = 3

                    if H < 0:
                        H = 0

                    H_id = round(H)
                    my_AOD_id = h2aod_id(H)

                    if my_Mu0 > 80:
                        my_Mu0 = 80

                    if my_Mu0 < 0:
                        my_Mu0 = 0

                    if my_Mu > 72:
                        my_Mu = 72

                    if my_Mu < 0:
                        my_Mu = 0

                    if my_RAA > 180:
                        my_RAA = 180

                    if my_RAA < 0:
                        my_RAA = 0

                    Pa_PList[i][j] = Pa[int(my_AOD_id)][int(H_id)][int(my_Mu0 / 2)][int(my_Mu / 2)][int(my_RAA / 5)]
                    Ra_PList[i][j] = Ra[int(my_AOD_id)][int(H_id)][int(my_Mu0 / 2)][int(my_Mu / 2)][int(my_RAA / 5)]
                    TT_PList[i][j] = TT[int(my_AOD_id)][int(H_id)][int(my_Mu0 / 2)][int(my_Mu / 2)][int(my_RAA / 5)]

                    Rs[i][j] = (data[i][j] - Pa_PList[i][j]) / (
                                TT_PList[i][j] + Ra_PList[i][j] * (data[i][j] - Pa_PList[i][j]))

                if data[i][j] == -999:
                    Rs[i][j] = -999

            if cloud[i][j] == 1:
                Rs[i][j] = -999

    return Rs


def vi(NIR_f, R_f):
    NIR = gdal.Open(NIR_f)
    R = gdal.Open(R_f)
    NIR = NIR.ReadAsArray()
    R = R.ReadAsArray()
    ndvi = zs.NDVI(NIR, R)

    band_tif = gdal.GetDriverByName("GTiff")
    band_filename = './dataset/data/preprocess/result/ndsi.tif'
    band_dataset = band_tif.Create(band_filename, ndvi.shape[1], ndvi.shape[0], 1, gdal.GDT_Float32)
    band0 = band_dataset.GetRasterBand(1)
    band0.WriteArray(ndvi)


# 将文本的后xml后缀转换为vrt
def renaming(file):
    ext = os.path.splitext(file)
    if ext[1] == '.xml':
        new_name = ext[0] + '.vrt'
        if os.path.exists(new_name):
            os.remove(new_name)
        os.rename(file, new_name)


# 高程与气溶胶之间转换
def h2aod_id(h):
    if h < 1:
        aod_id = 2
    if 1 <= h <= 2:
        aod_id = 1
    if h > 2:
        aod_id = 0
    return aod_id
