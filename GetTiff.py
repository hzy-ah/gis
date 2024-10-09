import matplotlib.pyplot as plt
import netCDF4 as nc
from osgeo import gdal
import numpy as np
import xml.etree.ElementTree as ET
from pyhdf.SD import *
import h5py
import read_data as rd
import create_data as cd

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
def GetDataTiff_n(data, outputPath):
    driver = gdal.GetDriverByName('GTiff')
    Data = driver.Create(outputPath + 'data.tif', data.shape[2], data.shape[1], data.shape[0], gdal.GDT_Float32)
    for m in range(data.shape[0]):
        data_band = Data.GetRasterBand(m + 1)
        data_band.WriteArray(data)
    del Data


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
    vrtFile = gdal.Translate(VRT_path + 'AQUA.vrt', data, format='vrt')


# 修改VRT文件
def modify_VRT(VRT_path):
    tree = ET.parse(VRT_path + 'AQUA.vrt')
    root = tree.getroot()
    Metadata = root.find('Metadata')
    for MDI in Metadata.findall('MDI'):
        Metadata.remove(MDI)
    Metadata.set('domain', "GEOLOCATION")
    metadata = {
        "SRS": 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]',
        'X_DATASET': VRT_path + 'longitude.tif',
        'X_BAND': '1',
        # 'PIXEL_OFFSET': '0',
        # 'PIXEL_STEP': '1',
        'PIXEL_OFFSET': '0',
        'PIXEL_STEP': '5',
        'Y_DATASET': VRT_path + 'latitude.tif',
        'Y_BAND': '1',
        # 'LINE_OFFSET': '0',
        # 'LINE_STEP': '1'
        'LINE_OFFSET': '0',
        'LINE_STEP': '5'
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
                            dstNodata=Nodata, srcNodata=Nodata,
                            geoloc=True, resampleAlg=gdal.GRIORA_Bilinear)
    data = correctData.ReadAsArray()
    correctData = None
    return data
# data = h5py.File('F:/data/npp/04/RNSCA-RVIRS_npp_d20230408_t1345_svi01.h5')
# allkey = []
# for key in data.keys():
#   allkey.append(key)
# allkey2 = []
# for key in data[allkey[0]].keys():
#   allkey2.append(key)
# Ref = data[allkey[0]][allkey2[0]]['Reflectance']
# LST = np.array(Ref).astype(float)
# LST = rd.viirs_i(LST)
# LST = np.array(nc.Dataset('F:/Merge/Data/NPP/TOA(NPP-VIIRS)_2023_04_10.nc').groups['TOA_I'].variables['Top_of_Atmosphere_Reflectance'][:])[1]
# # LST = np.array(nc.Dataset('F:/new/fy3d/TOA(Terra-MODIS)_2023_04_16_0.nc').groups['TOA_500'].variables['Top_of_Atmosphere_Reflectance_500'][:])[2]
# print(LST.shape)
# # for i in range(LST.shape[0]):
# #     for j in range(LST.shape[1]):
# #         if LST[i, j] == 65535:
# #             LST[i, j] = 0
# driver = gdal.GetDriverByName('GTiff')
# Data = driver.Create('F:/new/jpss_own/toa_npp2.tif', LST.shape[1], LST.shape[0], 1, gdal.GDT_Float32)
# LST_band = Data.GetRasterBand(1)
# LST_band.WriteArray(LST)

# dataset = gdal.Open('F:/GPP/GPP.tif')
# band = dataset.GetRasterBand(1)
# LST = band.ReadAsArray()
#
# # # # TOA
LST = np.array((nc.Dataset('F:/Merge/Data/Mon_FY3D/Mon-GPP(FY3D-MERSI)_2023_05.nc').groups['GPP'].variables['GPP'][:]/30) )
# # LST = np.array(nc.Dataset('F:/new/fy3d/gf/SRF(FY3D-MERSI)_2023_04_16.nc').groups['SRF_500'].variables['Surface_Reflectance_500'][:])[2]
# print(LST[0])
for i in range(LST.shape[0]):
    for j in range(LST.shape[1]):
        if LST[i][j] < 0:
            LST[i][j] = 0
        if np.isnan(LST[i][j]) == 1:
            LST[i][j] = 0



LST = cd.replace_fill_value(LST)
LST = cd.replace_fill_value(LST)
print(LST.shape)

driver = gdal.GetDriverByName('GTiff')
Data = driver.Create('F:/GPP/gpp_own.tif', LST.shape[1], LST.shape[0], 1, gdal.GDT_Float32)
LST_band = Data.GetRasterBand(1)
LST_band.WriteArray(LST)
# #
# # 画图
# plt.imshow(LST)
# plt.show()

# NDVI
# LST = np.array(nc.Dataset('F:/new/fy3d/gf/NDVI(FY3D-MERSI)_2023_04_16.nc').groups['NDVI'].variables['NDVI'][:])
# # LST = np.array(nc.Dataset('F:/Merge/Data/JPSS1/NDVI(JPSS1-VIIRS)_2023_05_22.nc').groups['NDVI'].variables['NDVI'][:])
# print(LST.shape)
# # for i in range(LST.shape[0]):
# #     for j in range(LST.shape[1]):
# #         if LST[i, j] == 65535:
# #             LST[i, j] = 0
# driver = gdal.GetDriverByName('GTiff')
# Data = driver.Create('F:/new/fy3d/gf/ndvi.tif', LST.shape[1], LST.shape[0], 1, gdal.GDT_Float32)
# LST_band = Data.GetRasterBand(1)
# LST_band.WriteArray(LST)
#
# # 画图
# plt.imshow(LST)
# plt.show()

# # LST
# LST = np.array(nc.Dataset('F:/Merge/Data/FY3D/LST(FY3D-MERSI)2023-04-16.nc').groups['LST'].variables['LST'][:])
# LST = np.array(nc.Dataset('F:/new/fy3d/LST(Terra-MODIS)_2023_04_16.nc').groups['LST'].variables['LST'][:])
# print(LST.shape)
# # for i in range(LST.shape[0]):
# #     for j in range(LST.shape[1]):
# #         if LST[i, j] == 65535:
# #             LST[i, j] = 0
# driver = gdal.GetDriverByName('GTiff')
# Data = driver.Create('f:/new/fy3d/LST_modis.tif', LST.shape[1], LST.shape[0], 1, gdal.GDT_Float32)
# LST_band = Data.GetRasterBand(1)
# LST_band.WriteArray(LST)
#
# # 画图
# plt.imshow(LST)
# plt.show()

# # BT
# LST = np.array(nc.Dataset('TOA(Terra-MODIS)_2023_04_16_0.nc').groups['BT_1000'].variables['Brightness_Temperature_1000'][:])[0]
# print(LST.shape)
# # for i in range(LST.shape[0]):
# #     for j in range(LST.shape[1]):
# #         if LST[i, j] == 65535:
# #             LST[i, j] = 0
# driver = gdal.GetDriverByName('GTiff')
# Data = driver.Create('./temp/BAND31_0.tif', LST.shape[1], LST.shape[0], 1, gdal.GDT_Float32)
# LST_band = Data.GetRasterBand(1)
# LST_band.WriteArray(LST)
#
# # 画图
# plt.imshow(LST)
# plt.show()

# # SRF
# LST = np.array(nc.Dataset('F:/new/09/LST(GK2A-AMI)_2023_05_26_05_50.nc').groups['LST'].variables['LST'][:])
# # LST = np.array(nc.Dataset('F:/Merge/Data/FY3D/04/LST(FY3D-MERSI)_2023_04_16.nc').groups['LST'].variables['LST'][:])
# # LST = np.array(nc.Dataset('F:/Merge/Data/JPSS1/SRF(JPSS1-VIIRS)_2023_05_22.nc').groups['SRF_I'].variables['Surface_Reflectance'][:])[0]
# print(LST.shape)
# # for i in range(LST.shape[0]):
# #     for j in range(LST.shape[1]):
# #         if LST[i, j] == 65535:
# #             LST[i, j] = 0
# driver = gdal.GetDriverByName('GTiff')
# Data = driver.Create('F:/new/09/lst.tif', LST.shape[1], LST.shape[0], 1, gdal.GDT_Float32)
# LST_band = Data.GetRasterBand(1)
# LST_band.WriteArray(LST)
#
# # 画图
# plt.imshow(LST)
# plt.show()


# isCloud = np.array(nc.Dataset('TOA(Aqua-MODIS)_2023_05_02.nc').groups['isCloud'].variables['isCloud'][:])
# print(isCloud.shape)
# SRF = np.array(nc.Dataset('SRF(Aqua-MODIS)_2023_05_02.nc').groups['SRF'].variables['Surface_Reflectance'][:])
# print(SRF.shape)
# for i in range(SRF.shape[0]):
#     for j in range(SRF.shape[1]):
#         for k in range(SRF.shape[2]):
#             if isCloud[0][j][k] == 1:
#                 SRF[i][j][k] = -2
# driver = gdal.GetDriverByName('GTiff')
# Data = driver.Create('5_2.tif', SRF.shape[2], SRF.shape[1], SRF.shape[0], gdal.GDT_Float32)
# for i in range(SRF.shape[0]):
#     data_band = Data.GetRasterBand(i + 1)
#     data_band.WriteArray(SRF[i])

#
# hdf = SD('./LST/MOD11_L2.A2023106.0205.061.2023107081507.hdf')
# datas0 = hdf.datasets()
# print(datas0)
# for i in datas0:
#     print(i)
# lon = hdf.select('Longitude')[:]
# lat = hdf.select('Latitude')[:]
# GetLocationTiff(lon, lat, './LST/')
# LST_B = np.array(hdf.select('LST')[:] * 0.02)
# # SRF_B = np.array(hdf.select('BAND32')[:] * 0.01)
# print(LST_B)
# GetDataTiff(LST_B, './LST/')
# Create_VRT('./LST/data.tif', './LST/')
# modify_VRT('./LST/')
# GeoCorrection('./LST/2.tif', './LST/AQUA.vrt', 0.01, 0)


# LST_1 = gdal.Open('./LST/1.tif').ReadAsArray()
# LST_2 = gdal.Open('./LST/2.tif').ReadAsArray()
# LST_B = np.zeros((LST_1.shape[0], LST_1.shape[1]))
# for i in range(LST_1.shape[0]):
#     for j in range(LST_1.shape[1]):
#         if LST_1[i][j] != 0 and LST_2[i][j] != 0:
#             LST_B[i][j] = (LST_1[i][j] + LST_2[i][j]) / 2
#         elif LST_1[i][j] == 0 and LST_2[i][j] != 0:
#             LST_B[i][j] = LST_2[i][j]
#         elif LST_1[i][j] != 0 and LST_2[i][j] == 0:
#             LST_B[i][j] = LST_1[i][j]
#         elif LST_1[i][j] == 0 and LST_2[i][j] == 0:
#             LST_B[i][j] = 0
# driver = gdal.GetDriverByName('GTiff')
# Data = driver.Create('./LST/LST_B.tif', LST_B.shape[1], LST_B.shape[0], 1, gdal.GDT_Float32)
# LST_band = Data.GetRasterBand(1)
# LST_band.WriteArray(LST_B)
#
# # 画图
# plt.imshow(LST_B)
# plt.show()
