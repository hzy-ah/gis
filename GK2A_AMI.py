#!/usr/bin/env python
############################
# GK2A L1B data Processing sample code2
#
# This program extracts user defined area's pixel value/latitude/
# longitude value from GK2A NetCDF4 file and converts digital count number to
# Albedo/Brightness Temperature.
# After that it saves converted data to new NetCDF4 file with geographic coordinates.
#
#
# Input  : GK2A L1B file [sample file: SW038/fd020ge] (netCDF4)
#		   GK2A conversion table(ASCII)
#
# process	: read input files -> cut user defined area from input
#			  -> convert digital count number to Albedo/Brightness Temperature
#			  -> save data with netCDF4 form
#
# Output : Albedo/Brightness Temperature from cut area (netCDF4)
#
# The output netCDF4 file includes next datas
# -user defined area's line & column size
# -user defined area's image_pixel_value
# -latitude & longitude of every pixel in cut area
# -user defined area's line/column number of left upper point in original GEOS image array (in global attribute)
# -user defined area's line/column number of right lower point in original GEOS image array(in global attribute)
#
############################
# Library import & Function define
############################
# 1. Library import
############################
import os
import numpy as np
from osgeo import gdal
import netCDF4 as nc
import xml.etree.ElementTree as ET


##########################
# 2. Define function : Full disc GEOS image(Line, Column)	->	(Latitude, Longitude)
##########################
def latlon_from_lincol_geos(Resolution, Line, Column):
    degtorad = 3.14159265358979 / 180.0
    if (Resolution == 0.5):
        COFF = 11000.5
        CFAC = 8.170135561335742e7
        LOFF = 11000.5
        LFAC = 8.170135561335742e7
    elif (Resolution == 1.0):
        COFF = 5500.5
        CFAC = 4.0850677806678705e7
        LOFF = 5500.5
        LFAC = 4.0850677806678705e7
    else:
        COFF = 2750.5
        CFAC = 2.0425338903339352e7
        LOFF = 2750.5
        LFAC = 2.0425338903339352e7
    sub_lon = 128.2
    sub_lon = sub_lon * degtorad

    x = degtorad * ((Column - COFF) * 2 ** 16 / CFAC)
    y = degtorad * ((Line - LOFF) * 2 ** 16 / LFAC)
    Sd = np.sqrt(
        np.abs((42164.0 * np.cos(x) * np.cos(y)) ** 2 - (np.cos(y) ** 2 + 1.006739501 * np.sin(y) ** 2) * 1737122264))
    Sn = (42164.0 * np.cos(x) * np.cos(y) - Sd) / (np.cos(y) ** 2 + 1.006739501 * np.sin(y) ** 2)
    S1 = 42164.0 - (Sn * np.cos(x) * np.cos(y))
    S2 = Sn * (np.sin(x) * np.cos(y))
    S3 = -Sn * np.sin(y)
    Sxy = np.sqrt(((S1 * S1) + (S2 * S2)))

    nlon = (np.arctan(S2 / S1) + sub_lon) / degtorad
    nlat = np.arctan((1.006739501 * S3) / Sxy) / degtorad

    return (nlat, nlon)


##########################
# 3. Define function : (Latitude, Longitude)	->	Full disc GEOS image(Line, Column)
##########################
def lincol_from_latlon_geos(Resolution, Latitude, Longitude):
    degtorad = 3.14159265358979 / 180.0
    if (Resolution == 0.5):
        COFF = 11000.5
        CFAC = 8.170135561335742e7
        LOFF = 11000.5
        LFAC = 8.170135561335742e7
    elif (Resolution == 1.0):
        COFF = 5500.5
        CFAC = 4.0850677806678705e7
        LOFF = 5500.5
        LFAC = 4.0850677806678705e7
    else:
        COFF = 2750.5
        CFAC = 2.0425338903339352e7
        LOFF = 2750.5
        LFAC = 2.0425338903339352e7

    sub_lon = 128.2
    sub_lon = sub_lon * degtorad
    Latitude = Latitude * degtorad
    Longitude = Longitude * degtorad

    c_lat = np.arctan(0.993305616 * np.tan(Latitude))
    RL = 6356.7523 / np.sqrt(1.0 - 0.00669438444 * np.cos(c_lat) ** 2.0)
    R1 = 42164.0 - RL * np.cos(c_lat) * np.cos(Longitude - sub_lon)
    R2 = -RL * np.cos(c_lat) * np.sin(Longitude - sub_lon)
    R3 = RL * np.sin(c_lat)
    Rn = np.sqrt(R1 ** 2.0 + R2 ** 2.0 + R3 ** 2.0)

    x = np.arctan(-R2 / R1) / degtorad
    y = np.arcsin(-R3 / Rn) / degtorad
    ncol = COFF + (x * 2.0 ** (-16) * CFAC)
    nlin = LOFF + (y * 2.0 ** (-16) * LFAC)
    return (nlin, ncol)


##########################
# 4. Define function : Cut image_pixel_values/lat/lon array with latitude, longitude from GEOS data array
#
# Input Argument
#  -Array: GEOS full disc image_pixel_values/latitude/longitude Array [array/numpy array]
#  -Resolution: GEOS data's Resolution(km) [float]
#  -Latitude1: Left upper position's latitude of user defined area (degree) [float]
#  -Longitude1: Left upper position's longitude of user defined area (degree) [float]
#  -Latitude2: Right lower position's latitude of user defined area (degree) [float]
#  -Longitude2: Right lower position's longitude of user defined area (degree) [float]
#
# Latitude1 >= Latitude2
# Longitude1 <= Latitude2
#
# Output: image_pixel_value/latitude/longitude array [numpy array]
##########################
def cut_with_latlon_geos(Array, Resolution, Latitude1, Longitude1, Latitude2, Longitude2, channel):
    Array = np.array(Array)
    if (Resolution == 0.5):
        Index_max = 22000
    elif (Resolution == 1.0):
        Index_max = 11000
    else:
        Index_max = 5500

    (Lin1, Col1) = lincol_from_latlon_geos(Resolution, Latitude1, Longitude1)
    (Lin2, Col2) = lincol_from_latlon_geos(Resolution, Latitude2, Longitude2)
    Col1 = int(np.floor(Col1))
    Lin1 = int(np.floor(Lin1))
    Col2 = int(np.ceil(Col2))
    Lin2 = int(np.ceil(Lin2))
    if channel == 'VI006':
        Col1 = Col1 - 8222
        Lin1 = Lin1 - 4696
        Col2 = Col2 - 8222
        Lin2 = Lin2 - 4696

    cut = np.zeros((Index_max, Index_max))
    if ((Col1 <= Col2) and (Lin1 <= Lin2) and (0 <= Col1) and (Col2 < Index_max) and (0 <= Lin1) and (
            Lin2 < Index_max)):
        cut = Array[Lin1:Lin2, Col1:Col2]

    return cut


def ReadData(input_ncfile_path, CT_path):
    ############################
    # Main Program Start
    ############################
    # 5. Input data path setup
    ############################

    left_upper_lat = 31
    left_upper_lon = 113
    right_lower_lat = 20
    right_lower_lon = 124

    ############################
    # 6. GK2A sample data file read
    # 样例文件读取
    ############################
    input_ncfile = nc.Dataset(input_ncfile_path, 'r', format='netcdf4')
    ipixel = input_ncfile.variables['image_pixel_values']
    channel = ipixel.getncattr('channel_name')
    # ipixel = input_ncfile.variables['image_pixel_values'][2348:3363, 4111:5067]

    ##########################
    # 7. Calculate latitude & longitude from GEOS image
    # 根据GEOS图像计算经纬度
    ##########################
    i = np.arange(0, input_ncfile.getncattr('number_of_columns'), dtype='f')
    j = np.arange(0, input_ncfile.getncattr('number_of_lines'), dtype='f')
    i, j = np.meshgrid(i, j)

    (geos_lat, geos_lon) = latlon_from_lincol_geos(1.0, j, i)

    ##########################
    # 8. Cut user defined area from GEOS image
    # 从GEOS图像中剪切用户定义的区域
    ##########################
    cut_pixel = cut_with_latlon_geos(ipixel[:], 1.0, left_upper_lat, left_upper_lon, right_lower_lat, right_lower_lon,
                                     channel)
    cut_lat = cut_with_latlon_geos(geos_lat, 1.0, left_upper_lat, left_upper_lon, right_lower_lat, right_lower_lon,
                                   channel)
    cut_lon = cut_with_latlon_geos(geos_lon, 1.0, left_upper_lat, left_upper_lon, right_lower_lat, right_lower_lon,
                                   channel)

    (ulc_lin, ulc_col) = lincol_from_latlon_geos(1.0, left_upper_lat, left_upper_lon)
    (lrc_lin, lrc_col) = lincol_from_latlon_geos(1.0, right_lower_lat, right_lower_lon)
    print(ulc_lin, ulc_col, lrc_lin, lrc_col)

    ############################
    # 9. image_pixel_values DQF processing
    # 异常值处理
    ############################
    cut_pixel[cut_pixel > 49151] = 0  # set error pixel's value to 0

    ############################
    # 10. image_pixel_values Bit Size per pixel masking
    # image_pixel_values每个像素掩码的比特深度
    ############################
    if ((channel == 'VI004') or (channel == 'VI005') or (channel == 'NR016')):
        mask = 0b0000011111111111  # 11bit mask
    elif ((channel == 'VI006') or (channel == 'NR013') or (channel == 'WV063')):
        mask = 0b0000111111111111  # 12bit mask
    elif (channel == 'SW038'):
        mask = 0b0011111111111111  # 14bit mask
    else:
        mask = 0b0001111111111111  # 13bit mask

    cut_pixel_masked = np.bitwise_and(cut_pixel, mask)
    print(channel)
    print(cut_pixel_masked)
    ############################
    # 11. image pixel value -> Albedo/Brightness Temperature
    # 像素值 -> 转化为反射率/亮温
    ############################
    AL_postfix = '_con_alb.txt'
    BT_postfix = '_con_bt.txt'
    if (channel[0:2] == 'VI') or (channel[0:2] == 'NR'):
        conversion_table = np.loadtxt(CT_path + channel + AL_postfix, 'float64')
        convert_data = 'albedo'
    else:
        conversion_table = np.loadtxt(CT_path + channel + BT_postfix, 'float64')
        convert_data = 'brightness_temperature'

    cut_pixel_masked_converted = conversion_table[cut_pixel_masked]  # pixel data : table value / 1:1 matching

    input_ncfile.close()
    print(cut_pixel_masked_converted)
    return cut_pixel_masked_converted, cut_lat, cut_lon, channel


# 获取VI006波段的数据
def ReadData_VI006(input_ncfile_path, CT_path):
    left_upper_lat = 31
    left_upper_lon = 113
    right_lower_lat = 20
    right_lower_lon = 124

    input_ncfile = nc.Dataset(input_ncfile_path, 'r', format='netcdf4')
    ipixel = input_ncfile.variables['image_pixel_values'][4696:6728, 8222:10136]
    channel = input_ncfile.variables['image_pixel_values'].getncattr('channel_name')
    # ipixel = input_ncfile.variables['image_pixel_values']

    ##########################
    # 7. Calculate latitude & longitude from GEOS image
    # 根据GEOS图像计算经纬度
    ##########################
    i = np.arange(0, input_ncfile.getncattr('number_of_columns'), dtype='f')
    j = np.arange(0, input_ncfile.getncattr('number_of_lines'), dtype='f')
    i, j = np.meshgrid(i, j)
    i = i[4696:6728, 8222:10136]
    j = j[4696:6728, 8222:10136]

    # i = np.arange(0, ipixel.shape[1], dtype='f')
    # j = np.arange(0, ipixel.shape[0], dtype='f')
    # i, j = np.meshgrid(i, j)

    (geos_lat, geos_lon) = latlon_from_lincol_geos(0.5, j, i)

    cut_pixel = cut_with_latlon_geos(ipixel, 0.5, left_upper_lat, left_upper_lon, right_lower_lat, right_lower_lon,
                                     channel)
    cut_lat = cut_with_latlon_geos(geos_lat, 0.5, left_upper_lat, left_upper_lon, right_lower_lat, right_lower_lon,
                                   channel)
    cut_lon = cut_with_latlon_geos(geos_lon, 0.5, left_upper_lat, left_upper_lon, right_lower_lat, right_lower_lon,
                                   channel)

    cut_pixel[cut_pixel > 49151] = 0

    mask = 0b0000111111111111

    cut_pixel_masked = np.bitwise_and(cut_pixel, mask)

    AL_postfix = '_con_alb.txt'
    BT_postfix = '_con_bt.txt'
    conversion_table = np.loadtxt(CT_path + channel + AL_postfix, 'float64')
    convert_data = 'albedo'

    cut_pixel_masked_converted = conversion_table[cut_pixel_masked]

    input_ncfile.close()
    return cut_pixel_masked_converted, cut_lat, cut_lon


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
    Data = driver.Create(outputPath + 'data.tif', data.shape[1], data.shape[0], 1, gdal.GDT_Float32)
    data_band = Data.GetRasterBand(1)
    data_band.WriteArray(data)
    del Data


# 生成VRT文件
def Create_VRT(dataset, VRT_path):
    data = gdal.Open(dataset)
    vrtFile = gdal.Translate(VRT_path + 'GK2A.vrt', data, format='vrt')


# 修改VRT文件
def modify_VRT(VRT_path):
    tree = ET.parse(VRT_path + 'GK2A.vrt')
    root = tree.getroot()
    Metadata = root.find('Metadata')
    for MDI in Metadata.findall('MDI'):
        Metadata.remove(MDI)
    Metadata.set('domain', "GEOLOCATION")
    metadata = {
        "SRS": 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AXIS["Lat",NORTH],AXIS["Long",EAST],AUTHORITY["EPSG","4326"]]',
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
    tree.write(VRT_path + 'GK2A.vrt')


# 几何校正
def GeoCorrection(TIF, VRT, Resolution, Nodata):
    correctData = gdal.Warp(TIF, VRT, dstSRS="EPSG:4490", format='GTiff',
                            xRes=Resolution, yRes=Resolution, outputBounds=(115, 22,
                                                                            122, 29),
                            dstNodata=Nodata, srcNodata=Nodata,
                            geoloc=True, resampleAlg=gdal.GRIORA_Bilinear)
    correctData = None


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


def GetData(fileName, filePath_TIF):
    CT_Path = './conversion_table/'
    data, latitude, longitude, channel = ReadData(fileName, CT_Path)
    GetLocationTiff(longitude, latitude, filePath_TIF)
    GetDataTiff(data, filePath_TIF)
    Create_VRT(filePath_TIF + 'data.tif', filePath_TIF)
    modify_VRT(filePath_TIF)
    output_file_path = filePath_TIF + '/Output/'
    if not os.path.exists(output_file_path):
        os.mkdir(output_file_path)

    GeoCorrection(output_file_path + 'GK2A_1000_' + channel + '.tif', filePath_TIF + 'GK2A.vrt', 0.01, 255)
    GeoCorrection(output_file_path + 'GK2A_2000_' + channel + '.tif', filePath_TIF + 'GK2A.vrt', 0.02, 255)
    TOA_1000 = gdal.Open(output_file_path + 'GK2A_1000_' + channel + '.tif').ReadAsArray()
    TOA_2000 = gdal.Open(output_file_path + 'GK2A_2000_' + channel + '.tif').ReadAsArray()
    lat_1000, lon_1000 = GetGeolocation(output_file_path + 'GK2A_1000_' + channel + '.tif')
    lat_2000, lon_2000 = GetGeolocation(output_file_path + 'GK2A_2000_' + channel + '.tif')
    return TOA_1000, TOA_2000, channel, lat_1000, lon_1000, lat_2000, lon_2000


def GetData_VI006(fileName, filePath_TIF):
    CT_Path = './conversion_table/'
    data, latitude, longitude = ReadData_VI006(fileName, CT_Path)
    GetLocationTiff(longitude, latitude, filePath_TIF)
    GetDataTiff(data, filePath_TIF)
    Create_VRT(filePath_TIF + 'data.tif', filePath_TIF)
    modify_VRT(filePath_TIF)
    output_file_path = filePath_TIF + '/Output/'
    GeoCorrection(output_file_path + 'GK2A_VI006_1000' + '.tif', filePath_TIF + 'GK2A.vrt', 0.01, 255)
    GeoCorrection(output_file_path + 'GK2A_VI006_2000' + '.tif', filePath_TIF + 'GK2A.vrt', 0.02, 255)
    TOA_1000 = gdal.Open(output_file_path + 'GK2A_VI006_1000' + '.tif').ReadAsArray()
    TOA_2000 = gdal.Open(output_file_path + 'GK2A_VI006_2000' + '.tif').ReadAsArray()
    return TOA_1000, TOA_2000


def GetAngle(filePath):
    NC = nc.Dataset(filePath)
    sun = NC.variables['sun_position']
    sc = NC.variables['sc_position']
    so_z = sun.getncattr('sun_zenith_angle')
    so_a = sun.getncattr('sun_azimuth_angle')
    sa_z = sc.getncattr('sc_zenith_angle')
    sa_a = sc.getncattr('sc_azimuth_angle')
    return so_z, so_a, sa_z, sa_a


def TOA_nc(TOA_1000, TOA_2000, BT_2000, lat_1000, lon_1000, lat_2000, lon_2000, isCloud_1000, isCloud_2000, out_dir,
           dateTime):
    out_path = out_dir + '/' + 'TOA(GK2A-AMI)_' + dateTime + '.nc'
    f_w = nc.Dataset(out_path, 'w', format='NETCDF4')
    f_w.createDimension('lat_1000', len(lat_1000))
    f_w.createDimension('lon_1000', len(lon_1000))
    f_w.createDimension('lat_2000', len(lat_2000))
    f_w.createDimension('lon_2000', len(lon_2000))
    f_w.createDimension('size1', TOA_1000.shape[0])
    f_w.createDimension('size2', TOA_2000.shape[0])
    f_w.createDimension('size3', BT_2000.shape[0])
    f_w.createDimension('size4', 1)
    f_w.createDimension('size5', 1)

    group1 = f_w.createGroup('TOA_1000')
    group1.createVariable('lat_1000', np.float32, ('lat_1000'))
    group1.variables['lat_1000'][:] = lat_1000
    group1.createVariable('lon_1000', np.float32, ('lon_1000'))
    group1.variables['lon_1000'][:] = lon_1000
    Top_of_Atmosphere_Reflectance_1000 = group1.createVariable('Top_of_Atmosphere_Reflectance_1000', np.float32,
                                                               ('size1', 'lat_1000', 'lon_1000'))
    group1.variables['Top_of_Atmosphere_Reflectance_1000'][:] = TOA_1000
    Top_of_Atmosphere_Reflectance_1000.datetime = dateTime
    Top_of_Atmosphere_Reflectance_1000.band_names = ['VI004', 'VI005', 'VI006', 'VI008']

    group2 = f_w.createGroup('TOA_2000')
    group2.createVariable('lat_2000', np.float32, ('lat_2000'))
    group2.variables['lat_2000'][:] = lat_2000
    group2.createVariable('lon_2000', np.float32, ('lon_2000'))
    group2.variables['lon_2000'][:] = lon_2000
    Top_of_Atmosphere_Reflectance_2000 = group2.createVariable('Top_of_Atmosphere_Reflectance_2000', np.float32,
                                                               ('size2', 'lat_2000', 'lon_2000'))
    group2.variables['Top_of_Atmosphere_Reflectance_2000'][:] = TOA_2000
    Top_of_Atmosphere_Reflectance_2000.datetime = dateTime
    Top_of_Atmosphere_Reflectance_2000.band_names = ['NR013', 'NR016']

    group3 = f_w.createGroup('BT_2000')
    group3.createVariable('lat_2000', np.float32, ('lat_2000'))
    group3.variables['lat_2000'][:] = lat_2000
    group3.createVariable('lon_2000', np.float32, ('lon_2000'))
    group3.variables['lon_2000'][:] = lon_2000
    Brightness_Temperature_2000 = group3.createVariable('Brightness_Temperature_2000', np.float32,
                                                        ('size3', 'lat_2000', 'lon_2000'))
    group3.variables['Brightness_Temperature_2000'][:] = BT_2000
    Brightness_Temperature_2000.datetime = dateTime
    Brightness_Temperature_2000.band_names = ['SW038', 'IR112', 'IR123']

    group4 = f_w.createGroup('isCloud_1000')
    group4.createVariable('lat_1000', np.float32, ('lat_1000'))
    group4.variables['lat_1000'][:] = lat_1000
    group4.createVariable('lon_1000', np.float32, ('lon_1000'))
    group4.variables['lon_1000'][:] = lon_1000
    Cloud = group4.createVariable('isCloud_1000', np.float32, ('size4', 'lat_1000', 'lon_1000'))
    group4.variables['isCloud_1000'][:] = isCloud_1000
    Cloud.datetime = dateTime
    Cloud.band_names = ['isCloud_1000']

    group5 = f_w.createGroup('isCloud_2000')
    group5.createVariable('lat_2000', np.float32, ('lat_2000'))
    group5.variables['lat_2000'][:] = lat_2000
    group5.createVariable('lon_2000', np.float32, ('lon_2000'))
    group5.variables['lon_2000'][:] = lon_2000
    Cloud = group5.createVariable('isCloud_2000', np.float32, ('size5', 'lat_2000', 'lon_2000'))
    group5.variables['isCloud_2000'][:] = isCloud_2000
    Cloud.datetime = dateTime
    Cloud.band_names = ['isCloud_2000']
