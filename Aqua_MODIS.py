import os
import numpy as np
from osgeo import gdal
from pyhdf.HDF import *
from pyhdf.SD import *
from pyhdf.VS import *
import netCDF4 as nc
import xml.etree.ElementTree as ET
import math


# 获取日间文件列表
def GetDaytimeFiles(file_list_inDate, L1B_Dir):
    file_list_day = list()
    for file in file_list_inDate:
        hdf = HDF(L1B_Dir + '/' + file)
        vs = hdf.vstart()
        Scans_Number = vs.attach('Number of Scans')[:][0][0]
        Scans_DayNumber = vs.attach('Number of Day mode scans')[:][0][0]
        Scans_NightNumber = vs.attach('Number of Night mode scans')[:][0][0]
        if Scans_Number == Scans_DayNumber:
            file_list_day.append(file)
    return file_list_day


# 生成经纬度TIFF文件
def GetLocationTiff(longitude, latitude, outputPath):
    driver = gdal.GetDriverByName('GTiff')
    # longitude
    lon = driver.Create(outputPath + 'longitude.tif', longitude.shape[1], longitude.shape[0], 1, gdal.GDT_Float32)
    lon_band = lon.GetRasterBand(1)
    lon_band.WriteArray(longitude)
    del lon
    # latitude
    lat = driver.Create(outputPath + 'latitude.tif', latitude.shape[1], latitude.shape[0], 1, gdal.GDT_Float32)
    lat_band = lat.GetRasterBand(1)
    lat_band.WriteArray(latitude)
    del lat


# 生成数据TIFF文件
def GetDataTiff(data, outputPath):
    driver = gdal.GetDriverByName('GTiff')
    Data = driver.Create(outputPath + 'data.tif', data.shape[2], data.shape[1], data.shape[0], gdal.GDT_Float32)
    for i in range(data.shape[0]):
        data_band = Data.GetRasterBand(i + 1)
        data_band.WriteArray(data[i])
    del Data


# 生成VRT文件
def Create_VRT(dataset, VRT_path):
    data = gdal.Open(dataset)
    vrtFile = gdal.Translate(VRT_path + 'AQUA.vrt', data, format='vrt')


# 修改VRT文件
def modify_VRT(VRT_path, mark):
    tree = ET.parse(VRT_path + 'AQUA.vrt')
    root = tree.getroot()
    Metadata = root.find('Metadata')
    for MDI in Metadata.findall('MDI'):
        Metadata.remove(MDI)
    Metadata.set('domain', "GEOLOCATION")
    if mark == True:
        metadata = {
            "SRS": 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]',
            'X_DATASET': VRT_path + 'longitude.tif',
            'X_BAND': '1',
            'PIXEL_OFFSET': '2',
            'PIXEL_STEP': '5',
            'Y_DATASET': VRT_path + 'latitude.tif',
            'Y_BAND': '1',
            'LINE_OFFSET': '2',
            'LINE_STEP': '5'
        }
    else:
        metadata = {
            "SRS": 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]',
            'X_DATASET': VRT_path + 'longitude.tif',
            'X_BAND': '1',
            'PIXEL_OFFSET': '0',
            'PIXEL_STEP': '1',
            'Y_DATASET': VRT_path + 'latitude.tif',
            'Y_BAND': '1',
            'LINE_OFFSET': '0',
            'LINE_STEP': '1'
        }
    for k in metadata.keys():
        MDI = ET.SubElement(Metadata, 'MDI')
        MDI.set('key', k)
        MDI.text = metadata[k]
    tree.write(VRT_path + 'AQUA.vrt')


# 几何校正
def GeoCorrection(TIF, VRT, Resolution, Nodata):
    correctData = gdal.Warp(TIF, VRT, dstSRS="EPSG:4490", format='GTiff',
                            xRes=Resolution, yRes=Resolution, outputBounds=(115, 22,
                                                                            122, 29),
                            srcNodata=Nodata,
                            geoloc=True, resampleAlg=gdal.GRIORA_Bilinear)
    data = correctData.ReadAsArray()
    correctData = None
    return data


# 获取NoData数据下标
def GetNodataFlag(fileName, filePath_TIF):
    hdf = SD(fileName)
    longitude0 = hdf.select('Longitude')[:]
    latitude0 = hdf.select('Latitude')[:]
    rows = longitude0.shape[0]
    cols = longitude0.shape[1]
    flag = []
    for i in range(rows):
        for j in range(cols):
            if longitude0[i][j] == -999 or latitude0[i][j] == -999:
                flag.append(i)
            break
    longitude = np.delete(longitude0, flag, axis=0)
    latitude = np.delete(latitude0, flag, axis=0)
    GetLocationTiff(longitude, latitude, filePath_TIF)
    return flag


# 单个文件经纬度范围判断
def isNeed(filePath_TIF):
    longitude = gdal.Open(filePath_TIF + 'longitude.tif').ReadAsArray()
    latitude = gdal.Open(filePath_TIF + 'latitude.tif').ReadAsArray()
    lat_min = min(np.hstack(latitude))
    lat_max = max(np.hstack(latitude))
    lon_min = min(np.hstack(longitude))
    lon_max = max(np.hstack(longitude))
    if lat_min > 29 or lat_max < 22 or lon_min > 122 or lon_max < 115:
        return False
    else:
        return True


# 获取辐射定标系数
def GetCalibrationCoefficient(hdf):
    calibration = np.zeros((10, 2))
    calibration[0][0] = hdf.select('EV_250_Aggr1km_RefSB').reflectance_scales[0]
    calibration[0][1] = hdf.select('EV_250_Aggr1km_RefSB').reflectance_offsets[0]
    calibration[1][0] = hdf.select('EV_250_Aggr1km_RefSB').reflectance_scales[1]
    calibration[1][1] = hdf.select('EV_250_Aggr1km_RefSB').reflectance_offsets[1]
    calibration[2][0] = hdf.select('EV_500_Aggr1km_RefSB').reflectance_scales[0]
    calibration[2][1] = hdf.select('EV_500_Aggr1km_RefSB').reflectance_offsets[0]
    calibration[3][0] = hdf.select('EV_500_Aggr1km_RefSB').reflectance_scales[1]
    calibration[3][1] = hdf.select('EV_500_Aggr1km_RefSB').reflectance_offsets[1]
    calibration[4][0] = hdf.select('EV_500_Aggr1km_RefSB').reflectance_scales[2]
    calibration[4][1] = hdf.select('EV_500_Aggr1km_RefSB').reflectance_offsets[2]
    calibration[5][0] = hdf.select('EV_500_Aggr1km_RefSB').reflectance_scales[3]
    calibration[5][1] = hdf.select('EV_500_Aggr1km_RefSB').reflectance_offsets[3]
    calibration[6][0] = hdf.select('EV_500_Aggr1km_RefSB').reflectance_scales[4]
    calibration[6][1] = hdf.select('EV_500_Aggr1km_RefSB').reflectance_offsets[4]
    calibration[7][0] = hdf.select('EV_1KM_RefSB').reflectance_scales[14]
    calibration[7][1] = hdf.select('EV_1KM_RefSB').reflectance_offsets[14]
    calibration[8][0] = hdf.select('EV_1KM_Emissive').radiance_scales[10]
    calibration[8][1] = hdf.select('EV_1KM_Emissive').radiance_offsets[10]
    calibration[9][0] = hdf.select('EV_1KM_Emissive').radiance_scales[11]
    calibration[9][1] = hdf.select('EV_1KM_Emissive').radiance_offsets[11]
    return calibration


# Aqua_MODIS辐射校正
def RadioCorrection(data1, data2, data3, calibration):
    tif1 = gdal.Open(data1)
    tif2 = gdal.Open(data2)
    tif3 = gdal.Open(data3)
    datasets1 = tif1.ReadAsArray()
    datasets2 = tif2.ReadAsArray()
    datasets3 = tif3.ReadAsArray()
    band_500 = np.zeros((datasets1.shape[0], datasets1.shape[1], datasets1.shape[2]))
    band_1000 = np.zeros((datasets2.shape[0], datasets2.shape[1], datasets2.shape[2]))
    band_1000_500 = np.zeros((datasets3.shape[0], datasets3.shape[1], datasets3.shape[2]))
    # band1,2,3,4,5,6,7
    for i in range(datasets1.shape[0]):
        for j in range(datasets1.shape[1]):
            for k in range(datasets1.shape[2]):
                if datasets1[i][j][k] != -999:
                    band_500[i][j][k] = calibration[i][0] * (datasets1[i][j][k] - calibration[i][1])
    # band26
    for j in range(datasets2.shape[1]):
        for k in range(datasets2.shape[2]):
            if datasets2[0][j][k] != -999:
                band_1000[0][j][k] = calibration[7][0] * (datasets2[0][j][k] - calibration[7][1])
    # band26_1000_500
    for j in range(datasets3.shape[1]):
        for k in range(datasets3.shape[2]):
            if datasets3[0][j][k] != -999:
                band_1000_500[0][j][k] = calibration[7][0] * (datasets3[0][j][k] - calibration[7][1])
    C1 = 1.19104356e+8
    C2 = 1.4387685e+4
    Cw31 = 11.03
    Cw32 = 12.02
    # band31
    for j in range(datasets2.shape[1]):
        for k in range(datasets2.shape[2]):
            if datasets2[1][j][k] != -999:
                I = calibration[8][0] * (datasets2[1][j][k] - calibration[8][1])
                band_1000[1][j][k] = C2 / (Cw31 * math.log(1 + C1 / (I * Cw31 ** 5)))
    # band31_1000_500
    for j in range(datasets3.shape[1]):
        for k in range(datasets3.shape[2]):
            if datasets3[1][j][k] != -999:
                I = calibration[8][0] * (datasets3[1][j][k] - calibration[8][1])
                band_1000_500[1][j][k] = C2 / (Cw31 * math.log(1 + C1 / (I * Cw31 ** 5)))
    # band32
    for j in range(datasets2.shape[1]):
        for k in range(datasets2.shape[2]):
            if datasets2[2][j][k] != -999:
                I = calibration[9][0] * (datasets2[2][j][k] - calibration[9][1])
                band_1000[2][j][k] = C2 / (Cw32 * math.log(1 + C1 / (I * Cw32 ** 5)))

    # band32_1000_500
    for j in range(datasets3.shape[1]):
        for k in range(datasets3.shape[2]):
            if datasets3[2][j][k] != -999:
                I = calibration[9][0] * (datasets3[2][j][k] - calibration[9][1])
                band_1000_500[2][j][k] = C2 / (Cw32 * math.log(1 + C1 / (I * Cw32 ** 5)))
    return band_500, band_1000, band_1000_500


# 单个文件经纬度获取
def GetGeolocation(data):
    datasets = gdal.Open(data)
    adfGeoTransform = datasets.GetGeoTransform()

    nXSize = int(datasets.RasterXSize)  # 列数
    nYSize = int(datasets.RasterYSize)  # 行数

    lat = np.ones((nYSize, nXSize))  # 用于存储每个像素的经度
    lon = np.ones((nYSize, nXSize))  # 用于存储每个像素的纬度
    for i in range(nYSize):
        for j in range(nXSize):
            px = adfGeoTransform[0] + i * adfGeoTransform[1] + j * adfGeoTransform[2]
            py = adfGeoTransform[3] + i * adfGeoTransform[4] + j * adfGeoTransform[5]
            lat[i][j] = py
            lon[i][j] = px
    lat_list = []
    lon_list = []
    for i in range(lat.shape[0]):
        lon_list.append(lon[i][0])
    for j in range(lon.shape[1]):
        lat_list.append(lat[0][j])
    return lat_list, lon_list


# 单个文件数据预处理
def GetData(fileName, flag, filePath_TIF, hour):
    hdf = SD(fileName)
    datas0 = hdf.datasets()
    flag_data = []
    for m in flag:
        for n in range(5):
            flag_data.append(m * 5 + n)
    getShape_geo = hdf.select('Longitude')[:]
    getShape_data = hdf.select('EV_1KM_RefSB')[:]
    data1 = np.ndarray(shape=[7, getShape_data.shape[1] - len(flag_data), getShape_data.shape[2]])
    data2 = np.ndarray(shape=[3, getShape_data.shape[1] - len(flag_data), getShape_data.shape[2]])
    angle = np.ndarray(shape=[4, getShape_geo.shape[0] - len(flag), getShape_geo.shape[1]])

    for channel in datas0:
        if channel == 'EV_1KM_RefSB':  # band_26
            data0 = hdf.select(channel)[:]
            data2[0] = np.delete(data0[14], flag_data, axis=0)
        elif channel == 'EV_1KM_Emissive':  # band_31,32
            data0 = hdf.select(channel)[:]
            data2[1] = np.delete(data0[10], flag_data, axis=0)
            data2[2] = np.delete(data0[11], flag_data, axis=0)
        elif channel == 'EV_250_Aggr1km_RefSB':  # band_1,2
            data0 = hdf.select(channel)[:]
            data1[0] = np.delete(data0[0], flag_data, axis=0)
            data1[1] = np.delete(data0[1], flag_data, axis=0)
        elif channel == 'EV_500_Aggr1km_RefSB':  # band_3,4,5,6,7
            data0 = hdf.select(channel)[:]
            data1[2] = np.delete(data0[0], flag_data, axis=0)
            data1[3] = np.delete(data0[1], flag_data, axis=0)
            data1[4] = np.delete(data0[2], flag_data, axis=0)
            data1[5] = np.delete(data0[3], flag_data, axis=0)
            data1[6] = np.delete(data0[4], flag_data, axis=0)
        elif channel == 'SensorZenith':
            data0 = hdf.select(channel)[:]
            angle[0] = np.delete(data0, flag, axis=0)
        elif channel == 'SensorAzimuth':
            data0 = hdf.select(channel)[:]
            angle[1] = np.delete(data0, flag, axis=0)
        elif channel == 'SolarZenith':
            data0 = hdf.select(channel)[:]
            angle[2] = np.delete(data0, flag, axis=0)
        elif channel == 'SolarAzimuth':
            data0 = hdf.select(channel)[:]
            angle[3] = np.delete(data0, flag, axis=0)

    # data1(TOA_500:1,2,3,4,5,6,7)
    GetDataTiff(data1, filePath_TIF)
    Create_VRT(filePath_TIF + 'data.tif', filePath_TIF)
    modify_VRT(filePath_TIF, True)
    output_file_path = filePath_TIF + '/Output/'
    if not os.path.exists(output_file_path):
        os.mkdir(output_file_path)
    GeoCorrection(output_file_path + 'Aqua_Data1_' + hour + '.tif', filePath_TIF + 'AQUA.vrt', 0.005, -999)
    # data2(TOA_1000(26),BT_1000(31,32)
    GetDataTiff(data2, filePath_TIF)
    Create_VRT(filePath_TIF + 'data.tif', filePath_TIF)
    modify_VRT(filePath_TIF, True)
    GeoCorrection(output_file_path + 'Aqua_Data2_' + hour + '.tif', filePath_TIF + 'AQUA.vrt', 0.01, -999)
    # data3((TOA_500(26),BT_500(31,32))
    GeoCorrection(output_file_path + 'Aqua_Data3_' + hour + '.tif', filePath_TIF + 'AQUA.vrt', 0.005, -999)
    # angle
    GetDataTiff(angle, filePath_TIF)
    Create_VRT(filePath_TIF + 'data.tif', filePath_TIF)
    modify_VRT(filePath_TIF, False)
    GeoCorrection(output_file_path + 'Aqua_angle_' + hour + '.tif', filePath_TIF + 'AQUA.vrt', 0.005, 25500)

    # 获取定标系数
    calibration = GetCalibrationCoefficient(hdf)
    # 辐射校正
    band_500, band_1000, band_1000_500 = RadioCorrection(output_file_path + 'Aqua_Data1_' + hour + '.tif',
                                                         output_file_path + 'Aqua_Data2_' + hour + '.tif',
                                                         output_file_path + 'Aqua_Data3_' + hour + '.tif',
                                                         calibration)
    # 获取经纬度
    lat_500, lon_500 = GetGeolocation(output_file_path + 'Aqua_Data1_' + hour + '.tif')
    lat_1000, lon_1000 = GetGeolocation(output_file_path + 'Aqua_Data2_' + hour + '.tif')
    angle_data = gdal.Open(output_file_path + 'Aqua_angle_' + hour + '.tif').ReadAsArray()
    sa_z = angle_data[0]
    sa_a = angle_data[1]
    so_z = angle_data[2]
    so_a = angle_data[3]
    return band_500, band_1000, band_1000_500, lat_500, lon_500, lat_1000, lon_1000, sa_z, sa_a, so_z, so_a


# 数据合并
def MergeData(band1, band2, NoData):
    band = band1
    for i in range(band1.shape[0]):
        for j in range(band2.shape[1]):
            if band1[i][j] == NoData and band2[i][j] == NoData:
                band[i][j] = NoData
            elif band1[i][j] == NoData and band2[i][j] != NoData:
                band[i][j] = band2[i][j]
            elif band1[i][j] != NoData and band2[i][j] == NoData:
                band[i][j] = band1[i][j]
            elif band1[i][j] != NoData and band2[i][j] != NoData:
                band[i][j] = (band1[i][j] + band2[i][j]) / 2
    return band


# 将TOA产品写入nc文件
def TOA_nc(band_500, band_1000, sa_a, sa_z, so_a, so_z, A, lat_500, lon_500, lat_1000, lon_1000, isCloud, Ins_Dir,
           date):
    Output_path = Ins_Dir + '/' + 'TOA(Aqua-MODIS)_' + date + '.nc'
    f_w = nc.Dataset(Output_path, 'w', format='NETCDF4')
    f_w.createDimension('lat_500', len(lat_500))
    f_w.createDimension('lon_500', len(lon_500))
    f_w.createDimension('lat_1000', len(lat_1000))
    f_w.createDimension('lon_1000', len(lon_1000))
    f_w.createDimension('size1', band_500.shape[0])
    f_w.createDimension('size2', 1)
    f_w.createDimension('size3', 2)

    # band1,2,3,4,5,6,7
    group1 = f_w.createGroup('TOA_500')

    group1.createVariable('lat_500', np.float32, ('lat_500'))
    group1.variables['lat_500'][:] = lat_500

    group1.createVariable('lon_500', np.float32, ('lon_500'))
    group1.variables['lon_500'][:] = lon_500

    TOA_500 = group1.createVariable('Top_of_Atmosphere_Reflectance_500', np.float32, ('size1', 'lat_500', 'lon_500'))
    group1.variables['Top_of_Atmosphere_Reflectance_500'][:] = band_500
    TOA_500.time = date
    TOA_500.band_names = [1, 2, 3, 4, 5, 6, 7]

    # band26
    group2 = f_w.createGroup('TOA_1000')

    group2.createVariable('lat_1000', np.float32, ('lat_1000'))
    group2.variables['lat_1000'][:] = lat_1000

    group2.createVariable('lon_1000', np.float32, ('lon_1000'))
    group2.variables['lon_1000'][:] = lon_1000

    TOA_1000 = group2.createVariable('Top_of_Atmosphere_Reflectance_1000', np.float32,
                                     ('size2', 'lat_1000', 'lon_1000'))
    group2.variables['Top_of_Atmosphere_Reflectance_1000'][:] = band_1000[0]
    TOA_1000.time = date
    TOA_1000.band_names = [26]

    # band31,32
    group3 = f_w.createGroup('BT_1000')

    group3.createVariable('lat_1000', np.float32, ('lat_1000'))
    group3.variables['lat_1000'][:] = lat_1000

    group3.createVariable('lon_1000', np.float32, ('lon_1000'))
    group3.variables['lon_1000'][:] = lon_1000

    BT_1000 = group3.createVariable('Brightness_Temperature_1000', np.float32,
                                    ('size3', 'lat_1000', 'lon_1000'))
    group3.variables['Brightness_Temperature_1000'][:] = band_1000[1:3]
    BT_1000.time = date
    BT_1000.band_names = [31, 32]
    # 云掩膜
    group4 = f_w.createGroup('isCloud')

    group4.createVariable('lat_500', np.float32, ('lat_500'))
    group4.variables['lat_500'][:] = lat_500

    group4.createVariable('lon_500', np.float32, ('lon_500'))
    group4.variables['lon_500'][:] = lon_500

    TOA_1000 = group4.createVariable('isCloud', np.float32,
                                     ('size2', 'lat_500', 'lon_500'))
    group4.variables['isCloud'][:] = isCloud
    TOA_1000.time = date
    TOA_1000.band_names = ['isCloud']

    # 角度
    group5 = f_w.createGroup('angle')

    group5.createVariable('lat_500', np.float32, ('lat_500'))
    group5.variables['lat_500'][:] = lat_500

    group5.createVariable('lon_500', np.float32, ('lon_500'))
    group5.variables['lon_500'][:] = lon_500

    soa = group5.createVariable('SolarAzimuth', np.float32,
                                ('size2', 'lat_500', 'lon_500'))
    group5.variables['SolarAzimuth'][:] = so_a
    soa.time = date
    soa.band_names = ['SolarAzimuth']

    soz = group5.createVariable('SolarZenith', np.float32,
                                ('size2', 'lat_500', 'lon_500'))
    group5.variables['SolarZenith'][:] = so_z
    soz.time = date
    soz.band_names = ['SolarZenith']

    saa = group5.createVariable('SensorAzimuth', np.float32,
                                ('size2', 'lat_500', 'lon_500'))
    group5.variables['SensorAzimuth'][:] = sa_a
    saa.time = date
    saa.band_names = ['SensorAzimuth']

    saz = group5.createVariable('SensorZenith', np.float32,
                                ('size2', 'lat_500', 'lon_500'))
    group5.variables['SensorZenith'][:] = sa_z
    saz.time = date
    saz.band_names = ['SensorZenith']

    a = group5.createVariable('RelativeAzimuth', np.float32,
                              ('size2', 'lat_500', 'lon_500'))
    group5.variables['RelativeAzimuth'][:] = A
    a.time = date
    a.band_names = ['RelativeAzimuth']

    f_w.close()


# 计算相对方位角
def relative_angle(sa_a, sa_z, so_a, so_z):
    # 角度转换参数
    dtr = math.pi / 180
    rtd = 180 / math.pi

    # 计算相对方位角
    A = np.zeros((sa_z.shape[0], sa_z.shape[1]))
    for i in range(sa_z.shape[0]):
        for j in range(sa_z.shape[1]):
            # A[i][j] = math.fabs(sa_a[i][j]-so_a[i][j])
            if sa_a[i][j] != 255 and so_a[i][j] != 255:
                raa = abs(sa_a[i][j] - so_a[i][j])
                if raa > 180:
                    raa = 360 - raa
                A[i][j] = math.acos(math.cos(sa_z[i][j] * dtr) * math.cos(so_z[i][j] * dtr) +
                                    math.sin(sa_z[i][j] * dtr) * math.sin(so_z[i][j] * dtr) *
                                    math.cos(raa * dtr)) * rtd
            else:
                A[i][j] = 180
    return A
