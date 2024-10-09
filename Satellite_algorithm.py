import re
import os
import shutil

import Aqua_MODIS
import Terra_MODIS
import NOAA19_AVHRR
import fmask_AVHRR
import GK2A_AMI
import Mon
import Year
import Mon_GPP_NPP
from fmask import *
from fmask_gk2a_night import *
from SRF import *
from LST import *
from Vegetation_Index import *


##########
# 生成数据TIFF文件
def Tiff(data, outputPath):
    driver = gdal.GetDriverByName('GTiff')
    Data = driver.Create(outputPath, data.shape[1], data.shape[0], 1, gdal.GDT_Float32)
    data_band = Data.GetRasterBand(1)
    data_band.WriteArray(data)
    del Data


# 生成数据TIFF文件
def DataTiff(data, outputPath):
    driver = gdal.GetDriverByName('GTiff')
    Data = driver.Create(outputPath, data.shape[2], data.shape[1], data.shape[0], gdal.GDT_Float32)
    for i in range(data.shape[0]):
        data_band = Data.GetRasterBand(i + 1)
        data_band.WriteArray(data[i])
    del Data


###########################

def Aqua_MODIS_InsProducts(L1B_Dir, Date, Ins_Dir):
    # 判断时间格式，不符合参数格式则程序退出！
    res = re.search('^[0-9]{4}-[0-9]{2}-[0-9]{2}$', Date)
    if res is None:
        print('Please check your input date!')
        return False
    # 获取年、月、日
    year = Date[0:4]
    month = Date[5:7]
    day = Date[8:10]
    # 获取所需日期的数据文件目录
    file_list_all = os.listdir(L1B_Dir)
    file_list_inDate = list()
    for file in file_list_all:
        if file.find(year + '_' + month + '_' + day) != -1:
            file_list_inDate.append(file)
    # 如果没有所需日期的数据文件，报错退出！
    if not file_list_inDate:
        print('Can not find the idea Aqua-MODIS L1B files at your input day!')
        return False
    # 去除非MODIS021KM的文件
    file_list_inDate_MOD021km = []
    for file in file_list_inDate:
        if file.find('MOD021KM') != -1:
            file_list_inDate_MOD021km.append(file)
    # 获取日间观测数据
    file_list_day_MOD021km = Aqua_MODIS.GetDaytimeFiles(file_list_inDate_MOD021km, L1B_Dir)
    # 若无日间观测数据，报错退出！
    if not file_list_day_MOD021km:
        print('Can not find the idea Aqua-MODIS daytime L1B Files!')
        return False
    # 检测是否存在temp文件夹，不存在则创建
    temp_file_path = Ins_Dir + '/temp/'
    if not os.path.exists(temp_file_path):
        os.mkdir(temp_file_path)
    # 处理日间数据文件->判断是否在经纬度范围内
    file_list_need = []
    for file_name in file_list_day_MOD021km:
        file_path = L1B_Dir + '/' + file_name
        flag = Aqua_MODIS.GetNodataFlag(file_path, temp_file_path)
        if Aqua_MODIS.isNeed(temp_file_path):
            file_list_need.append(file_name)
    # 处理日间数据文件->数据处理与合并
    file_num = len(file_list_need)
    band_1 = -999 * np.ones((1400, 1400))
    band_2 = -999 * np.ones((1400, 1400))
    band_3 = -999 * np.ones((1400, 1400))
    band_4 = -999 * np.ones((1400, 1400))
    band_5 = -999 * np.ones((1400, 1400))
    band_6 = -999 * np.ones((1400, 1400))
    band_7 = -999 * np.ones((1400, 1400))
    band_26_500 = -999 * np.ones((1400, 1400))
    band_31_500 = -999 * np.ones((1400, 1400))
    band_32_500 = -999 * np.ones((1400, 1400))
    band_26 = -999 * np.ones((700, 700))
    band_31 = -999 * np.ones((700, 700))
    band_32 = -999 * np.ones((700, 700))
    sa_z = 25500 * np.ones((1400, 1400))
    sa_a = 25500 * np.ones((1400, 1400))
    so_z = 25500 * np.ones((1400, 1400))
    so_a = 25500 * np.ones((1400, 1400))
    lat_500 = []
    lon_500 = []
    lat_1000 = []
    lon_1000 = []
    TOA_500 = np.ndarray(shape=[7, 1400, 1400])
    TOA_1000 = np.ndarray(shape=[3, 700, 700])
    TOA_1000_500 = np.ndarray(shape=[3, 1400, 1400])
    for file_name in file_list_need:
        print('数据预处理！')
        print(file_name)
        hour = file_name.split(year + '_' + month + '_' + day)[1][1:3]
        file_path = L1B_Dir + '/' + file_name
        flag = Aqua_MODIS.GetNodataFlag(file_path, temp_file_path)
        band_500, band_1000, band_1000_500, lat_500, lon_500, lat_1000, lon_1000, SA_Z, SA_A, SO_Z, SO_A = Aqua_MODIS.GetData(
            file_path, flag,
            temp_file_path, hour)
        band_1 = Aqua_MODIS.MergeData(band_1, band_500[0], -999)
        band_2 = Aqua_MODIS.MergeData(band_2, band_500[1], -999)
        band_3 = Aqua_MODIS.MergeData(band_3, band_500[2], -999)
        band_4 = Aqua_MODIS.MergeData(band_4, band_500[3], -999)
        band_5 = Aqua_MODIS.MergeData(band_5, band_500[4], -999)
        band_6 = Aqua_MODIS.MergeData(band_6, band_500[5], -999)
        band_7 = Aqua_MODIS.MergeData(band_7, band_500[6], -999)
        band_26 = Aqua_MODIS.MergeData(band_26, band_1000[0], -999)
        band_31 = Aqua_MODIS.MergeData(band_31, band_1000[1], -999)
        band_32 = Aqua_MODIS.MergeData(band_32, band_1000[2], -999)
        band_26_500 = Aqua_MODIS.MergeData(band_26_500, band_1000_500[0], -999)
        band_31_500 = Aqua_MODIS.MergeData(band_31_500, band_1000_500[1], -999)
        band_32_500 = Aqua_MODIS.MergeData(band_32_500, band_1000_500[2], -999)
        sa_z = Aqua_MODIS.MergeData(sa_z, SA_Z, 25500)
        sa_a = Aqua_MODIS.MergeData(sa_a, SA_A, 25500)
        so_z = Aqua_MODIS.MergeData(so_z, SO_Z, 25500)
        so_a = Aqua_MODIS.MergeData(so_a, SO_A, 25500)
        sa_z = sa_z / 100
        sa_a = sa_a / 100
        so_z = so_z / 100
        so_a = so_a / 100

    TOA_500[0] = band_1
    TOA_500[1] = band_2
    TOA_500[2] = band_3
    TOA_500[3] = band_4
    TOA_500[4] = band_5
    TOA_500[5] = band_6
    TOA_500[6] = band_7
    TOA_1000[0] = band_26
    TOA_1000[1] = band_31
    TOA_1000[2] = band_32
    TOA_1000_500[0] = band_26_500
    TOA_1000_500[1] = band_31_500
    TOA_1000_500[2] = band_32_500

    date = year + '_' + month + '_' + day
    # 删除文件夹
    # shutil.rmtree(temp_file_path)

    #############################################################################
    # 云检测临时输出
    # Tiff(band_26_500, temp_file_path + 'band_26_500.tif')
    # Tiff(band_31_500, temp_file_path + 'band_31_500.tif')
    # DataTiff(TOA_500, temp_file_path + 'TOA_500.tif')
    #############################################################################
    print('云掩膜！')
    # 云掩膜
    isCloud = Get_isCloud(TOA_500, band_26_500, band_31_500)
    print('相对方位角')
    A = Aqua_MODIS.relative_angle(sa_a, sa_z, so_a, so_z)
    print('TOA！')
    # TOA数据写入nc文件
    name = '(Aqua-MODIS)_'
    print('大气校正！')
    # 大气校正
    TOA_500_SRF = algorithm_SRF(TOA_500, sa_z, so_z, A, isCloud, lat_500, lon_500, Ins_Dir, name, date)
    print('植被指数')
    # 植被指数
    # RVI
    rvi = RVI(TOA_500_SRF[1], TOA_500_SRF[0])
    # NDVI
    ndvi = NDVI(TOA_500_SRF[1], TOA_500_SRF[0])
    # SAVI
    savi = SAVI(TOA_500_SRF[1], TOA_500_SRF[0], 0.5)
    # GVI
    gvi = GVI(TOA_500_SRF[2], TOA_500_SRF[0], TOA_500_SRF[5], TOA_500_SRF[6])
    # EVI
    evi = EVI(TOA_500_SRF[1], TOA_500_SRF[0], TOA_500_SRF[2])
    # LAI
    landcover = gdal.Open('./ndata/IGBP_lai_500.tif').ReadAsArray()
    lai = LAI(ndvi, landcover)

    # 云掩膜
    for i in range(isCloud.shape[0]):
        for j in range(isCloud.shape[1]):
            if isCloud[i][j] == 1:
                rvi[i][j] = -999
                ndvi[i][j] = -999
                savi[i][j] = -999
                gvi[i][j] = -999
                evi[i][j] = -999
                lai[i][j] = -999
    # 写入
    RVI_nc(rvi, lat_500, lon_500, Ins_Dir, name + date)
    NDVI_nc(ndvi, lat_500, lon_500, Ins_Dir, name + date)
    SAVI_nc(savi, lat_500, lon_500, Ins_Dir, name + date)
    GVI_nc(gvi, lat_500, lon_500, Ins_Dir, name + date)
    EVI_nc(evi, lat_500, lon_500, Ins_Dir, name + date)
    LAI_nc(lai, lat_500, lon_500, Ins_Dir, name + date)


def Terra_MODIS_InsProducts(L1B_Dir, Date, Ins_Dir):
    # 判断时间格式，不符合参数格式则程序退出！
    res = re.search('^[0-9]{4}-[0-9]{2}-[0-9]{2}$', Date)
    if res is None:
        print('Please check your input date!')
        return False
    # 获取年、月、日
    year = Date[0:4]
    month = Date[5:7]
    day = Date[8:10]
    # 获取所需日期的数据文件目录
    file_list_all = os.listdir(L1B_Dir)
    file_list_inDate = list()
    for file in file_list_all:
        if file.find(year + '_' + month + '_' + day) != -1:
            file_list_inDate.append(file)
    # 如果没有所需日期的数据文件，报错退出！
    if not file_list_inDate:
        print('Can not find the idea Terra-MODIS L1B files at your input day!')
        return False
    # 去除非MODIS021KM的文件
    file_list_inDate_MOD021km = []
    for file in file_list_inDate:
        if file.find('MOD021KM') != -1:
            file_list_inDate_MOD021km.append(file)
    # 获取日间观测数据
    file_list_day_MOD021km = Terra_MODIS.GetDaytimeFiles(file_list_inDate_MOD021km, L1B_Dir)
    # 若无日间观测数据，报错退出！
    if not file_list_day_MOD021km:
        print('Can not find the idea Terra-MODIS daytime L1B Files!')
        return False

    # 检测是否存在temp文件夹，不存在则创建
    temp_file_path = Ins_Dir + '/temp/'
    if not os.path.exists(temp_file_path):
        os.mkdir(temp_file_path)
    # 处理日间数据文件->判断是否在经纬度范围内
    file_list_need = []
    for file_name in file_list_day_MOD021km:
        file_path = L1B_Dir + '/' + file_name
        flag = Terra_MODIS.GetNodataFlag(file_path, temp_file_path)
        if Terra_MODIS.isNeed(temp_file_path):
            file_list_need.append(file_name)
    # 处理日间数据文件->数据处理与合并
    band_1 = -999 * np.ones((1400, 1400))
    band_2 = -999 * np.ones((1400, 1400))
    band_3 = -999 * np.ones((1400, 1400))
    band_4 = -999 * np.ones((1400, 1400))
    band_5 = -999 * np.ones((1400, 1400))
    band_6 = -999 * np.ones((1400, 1400))
    band_7 = -999 * np.ones((1400, 1400))
    band_26_500 = -999 * np.ones((1400, 1400))
    band_31_500 = -999 * np.ones((1400, 1400))
    band_32_500 = -999 * np.ones((1400, 1400))
    band_26 = -999 * np.ones((700, 700))
    band_31 = -999 * np.ones((700, 700))
    band_32 = -999 * np.ones((700, 700))
    sa_z = 25500 * np.ones((1400, 1400))
    sa_a = 25500 * np.ones((1400, 1400))
    so_z = 25500 * np.ones((1400, 1400))
    so_a = 25500 * np.ones((1400, 1400))
    lat_500 = []
    lon_500 = []
    lat_1000 = []
    lon_1000 = []
    TOA_500 = np.ndarray(shape=[7, 1400, 1400])
    TOA_1000 = np.ndarray(shape=[3, 700, 700])
    TOA_1000_500 = np.ndarray(shape=[3, 1400, 1400])
    for file_name in file_list_need:
        print('数据预处理！')
        print(file_name)
        hour = file_name.split(year + '_' + month + '_' + day)[1][1:3]
        file_path = L1B_Dir + '/' + file_name
        flag = Terra_MODIS.GetNodataFlag(file_path, temp_file_path)
        band_500, band_1000, band_1000_500, lat_500, lon_500, lat_1000, lon_1000, SA_Z, SA_A, SO_Z, SO_A = Terra_MODIS.GetData(
            file_path, flag,
            temp_file_path, hour)
        band_1 = Terra_MODIS.MergeData(band_1, band_500[0], -999)
        band_2 = Terra_MODIS.MergeData(band_2, band_500[1], -999)
        band_3 = Terra_MODIS.MergeData(band_3, band_500[2], -999)
        band_4 = Terra_MODIS.MergeData(band_4, band_500[3], -999)
        band_5 = Terra_MODIS.MergeData(band_5, band_500[4], -999)
        band_6 = Terra_MODIS.MergeData(band_6, band_500[5], -999)
        band_7 = Terra_MODIS.MergeData(band_7, band_500[6], -999)
        band_26 = Terra_MODIS.MergeData(band_26, band_1000[0], -999)
        band_31 = Terra_MODIS.MergeData(band_31, band_1000[1], -999)
        band_32 = Terra_MODIS.MergeData(band_32, band_1000[2], -999)
        band_26_500 = Terra_MODIS.MergeData(band_26_500, band_1000_500[0], -999)
        band_31_500 = Terra_MODIS.MergeData(band_31_500, band_1000_500[1], -999)
        band_32_500 = Terra_MODIS.MergeData(band_32_500, band_1000_500[2], -999)
        sa_z = Terra_MODIS.MergeData(sa_z, SA_Z, 25500)
        sa_a = Terra_MODIS.MergeData(sa_a, SA_A, 25500)
        so_z = Terra_MODIS.MergeData(so_z, SO_Z, 25500)
        so_a = Terra_MODIS.MergeData(so_a, SO_A, 25500)
        sa_z = sa_z / 100
        sa_a = sa_a / 100
        so_z = so_z / 100
        so_a = so_a / 100

    TOA_500[0] = band_1
    TOA_500[1] = band_2
    TOA_500[2] = band_3
    TOA_500[3] = band_4
    TOA_500[4] = band_5
    TOA_500[5] = band_6
    TOA_500[6] = band_7
    TOA_1000[0] = band_26
    TOA_1000[1] = band_31
    TOA_1000[2] = band_32
    TOA_1000_500[0] = band_26_500
    TOA_1000_500[1] = band_31_500
    TOA_1000_500[2] = band_32_500

    date = year + '_' + month + '_' + day
    # 删除文件夹
    # shutil.rmtree(temp_file_path)
    print('云掩膜！')
    # 云掩膜
    ######################################################
    # 云检测临时输出
    Tiff(band_26_500, temp_file_path + 'band_26_500.tif')
    Tiff(band_31_500, temp_file_path + 'band_31_500.tif')
    DataTiff(TOA_500, temp_file_path + 'TOA_500.tif')
    ####################################################
    isCloud = Get_isCloud(TOA_500, band_26_500, band_31_500)
    print('相对方位角')
    A = Terra_MODIS.relative_angle(sa_a, sa_z, so_a, so_z)
    print('TOA！')
    # TOA数据写入nc文件
    name = '(Terra-MODIS)_'
    print('大气校正！')
    # 大气校正
    TOA_500_SRF = algorithm_SRF(TOA_500, sa_z, so_z, A, isCloud, lat_500, lon_500, Ins_Dir, name, date)
    print('植被指数')
    # 植被指数
    # RVI
    rvi = RVI(TOA_500_SRF[1], TOA_500_SRF[0])
    # NDVI
    ndvi = NDVI(TOA_500_SRF[1], TOA_500_SRF[0])
    # SAVI
    savi = SAVI(TOA_500_SRF[1], TOA_500_SRF[0], 0.5)
    # GVI
    gvi = GVI(TOA_500_SRF[2], TOA_500_SRF[0], TOA_500_SRF[5], TOA_500_SRF[6])
    # EVI
    evi = EVI(TOA_500_SRF[1], TOA_500_SRF[0], TOA_500_SRF[2])
    # LAI
    landcover = gdal.Open('./ndata/IGBP_lai_500.tif').ReadAsArray()
    lai = LAI(ndvi, landcover)

    # 云掩膜
    for i in range(isCloud.shape[0]):
        for j in range(isCloud.shape[1]):
            if isCloud[i][j] == 1:
                rvi[i][j] = -999
                ndvi[i][j] = -999
                savi[i][j] = -999
                gvi[i][j] = -999
                evi[i][j] = -999
                lai[i][j] = -999
    # 写入
    RVI_nc(rvi, lat_500, lon_500, Ins_Dir, name + date)
    NDVI_nc(ndvi, lat_500, lon_500, Ins_Dir, name + date)
    SAVI_nc(savi, lat_500, lon_500, Ins_Dir, name + date)
    GVI_nc(gvi, lat_500, lon_500, Ins_Dir, name + date)
    EVI_nc(evi, lat_500, lon_500, Ins_Dir, name + date)
    LAI_nc(lai, lat_500, lon_500, Ins_Dir, name + date)


def NOAA19_AVHRR_InsProducts(L1B_Dir, Date, Ins_Dir):
    # 判断时间格式，不符合参数格式则程序退出！
    res = re.search('^[0-9]{4}-[0-9]{2}-[0-9]{2}$', Date)
    if res is None:
        print('Please check your input date!')
        return False
    # 获取年、月、日
    year = Date[0:4]
    month = Date[5:7]
    day = Date[8:10]
    # 获取所需日期的数据文件目录
    file_list_all = os.listdir(L1B_Dir)
    file_list_inDate = list()
    for file in file_list_all:
        if file.find(year + '_' + month + '_' + day) != -1:
            file_list_inDate.append(file)
    # 如果没有所需日期的数据文件，报错退出！
    if not file_list_inDate:
        print('Can not find the idea NOAA19-AVHRR L1B files at your input day!')
        return False
    # 获取日间观测数据
    file_list_day = NOAA19_AVHRR.GetDaytimeFiles(file_list_inDate, year + '_' + month + '_' + day)
    # 若无日间观测数据，报错退出！
    if not file_list_day:
        print('Can not find the idea NOAA19-AVHRR daytime L1B Files!')
        return False
    date = year + '_' + month + '_' + day
    # 检测是否存在temp文件夹，不存在则创建
    temp_file_path = Ins_Dir + '/temp/'
    if not os.path.exists(temp_file_path):
        os.mkdir(temp_file_path)
    band_1 = -999 * np.ones((700, 700))
    band_2 = -999 * np.ones((700, 700))
    band_3 = -999 * np.ones((700, 700))
    band_4 = -999 * np.ones((700, 700))
    band_5 = -999 * np.ones((700, 700))
    SA_Z = 25500 * np.ones((700, 700))
    SO_Z = 25500 * np.ones((700, 700))
    AAA = 25500 * np.ones((700, 700))
    lat = []
    lon = []
    TOA = np.ndarray(shape=[3, 700, 700])
    BT = np.ndarray(shape=[2, 700, 700])
    # 处理日间数据文件
    for file_name in file_list_day:
        print(file_name)
        print('数据预处理！')
        hour = file_name.split(year + '_' + month + '_' + day)[1][1:3]
        file_path = L1B_Dir + '/' + file_name
        Ref, Emissive, sa_z, so_z, A, latitude, longitude = NOAA19_AVHRR.noaa_read(file_path)
        if not NOAA19_AVHRR.isNeed(latitude, longitude):
            continue
        REF, EMISSIVE, lat, lon, sa_Z, so_Z, AA = NOAA19_AVHRR.GetData(Ref, Emissive, sa_z, so_z, A, latitude,
                                                                       longitude,
                                                                       temp_file_path, hour)

        band_1 = NOAA19_AVHRR.MergeData(band_1, REF[0], -999)
        band_2 = NOAA19_AVHRR.MergeData(band_2, REF[1], -999)
        band_3 = NOAA19_AVHRR.MergeData(band_3, REF[2], -999)
        band_4 = NOAA19_AVHRR.MergeData(band_4, EMISSIVE[0], -999)
        band_5 = NOAA19_AVHRR.MergeData(band_5, EMISSIVE[1], -999)
        SA_Z = NOAA19_AVHRR.MergeData(SA_Z, sa_Z, 25500)
        SO_Z = NOAA19_AVHRR.MergeData(SO_Z, so_Z, 25500)
        AAA = NOAA19_AVHRR.MergeData(AAA, AA, 25500)

    TOA[0] = band_1
    TOA[1] = band_2
    TOA[2] = band_3
    BT[0] = band_4
    BT[1] = band_5

    # 删除文件夹
    # shutil.rmtree(temp_file_path)
    # 云检测临时输出
    DataTiff(TOA, temp_file_path + 'TOA.tif')
    DataTiff(BT, temp_file_path + 'BT.tif')
    print('云掩膜！')
    # 云掩膜
    isCloud = fmask_AVHRR.Get_isCloud(TOA, BT)
    print('TOA！')
    # TOA数据写入nc文件
    name = '(NOAA19-AVHRR)_'
    print('大气校正！')
    # 大气校正
    TOA_SRF = algorithm_SRF_NOAA19(TOA, SA_Z, SO_Z, AAA, isCloud, lat, lon, Ins_Dir, name, date)
    print('植被指数')
    # 植被指数
    # RVI
    rvi = RVI(TOA_SRF[1], TOA_SRF[0])
    # NDVI
    ndvi = NDVI(TOA_SRF[1], TOA_SRF[0])
    # SAVI
    savi = SAVI(TOA_SRF[1], TOA_SRF[0], 0.5)
    # EVI
    evi = EVI(TOA_SRF[1], TOA_SRF[0], TOA_SRF[0])
    # LAI
    landcover = gdal.Open('./ndata/IGBP_lai_1000.tif').ReadAsArray()
    lai = LAI(ndvi, landcover)

    # 云掩膜
    for i in range(isCloud.shape[0]):
        for j in range(isCloud.shape[1]):
            if isCloud[i][j] == 1:
                rvi[i][j] = -999
                ndvi[i][j] = -999
                savi[i][j] = -999
                evi[i][j] = -999
                lai[i][j] = -999

    # 写入
    EVI_nc(evi, lat, lon, Ins_Dir, name + date)
    SAVI_nc(savi, lat, lon, Ins_Dir, name + date)
    NDVI_nc(ndvi, lat, lon, Ins_Dir, name + date)
    RVI_nc(rvi, lat, lon, Ins_Dir, name + date)
    LAI_nc(lai, lat, lon, Ins_Dir, name + date)


def GK2A_AMI_InsProducts(L1B_Dir, Date, Ins_Dir):
    # 判断时间格式，不符合参数格式则程序退出！
    res = re.search('^[0-9]{4}-[0-9]{2}-[0-9]{2}\s[0-9]{2}:[0-9]{2}$', Date)
    if res is None:
        print('Please check your input date and time!')
        return False
    # 获取年、月、日、小时、分钟
    year = Date[0:4]
    month = Date[5:7]
    day = Date[8:10]
    hour = Date[11:13]
    minute = Date[14:16]
    name = '(GK2A-AMI)_'
    dateTime = year + '_' + month + '_' + day + '_' + hour + '_' + minute
    # 获取所需日期时刻的数据文件目录
    file_list_all = os.listdir(L1B_Dir)
    file_list_inDate = list()
    for file in file_list_all:
        if file.find(year + month + day + hour + minute) != -1:
            file_list_inDate.append(file)
    # 如果没有所需日期时刻的数据文件，报错退出！
    if not file_list_inDate:
        print('Can not find the idea GK2A-AMI L1B file at your input date and time!')
        return False
    file_list_need = list()
    file_VI004 = ''
    file_VI005 = ''
    file_VI006 = ''
    file_VI008 = ''
    file_NR013 = ''
    file_NR016 = ''
    file_SW038 = ''
    file_IR112 = ''
    file_IR123 = ''
    for file in file_list_inDate:
        if file.find('vi004') != -1:
            file_VI004 = file
            file_list_need.append(file)
        elif file.find('vi005') != -1:
            file_VI005 = file
            file_list_need.append(file)
        elif file.find('vi008') != -1:
            file_VI008 = file
            file_list_need.append(file)
        elif file.find('nr013') != -1:
            file_NR013 = file
            file_list_need.append(file)
        elif file.find('nr016') != -1:
            file_NR016 = file
            file_list_need.append(file)
        elif file.find('sw038') != -1:
            file_SW038 = file
            file_list_need.append(file)
        elif file.find('ir112') != -1:
            file_IR112 = file
            file_list_need.append(file)
        elif file.find('ir123') != -1:
            file_IR123 = file
            file_list_need.append(file)
        elif file.find('vi006') != -1:
            file_VI006 = file
    # 判断是否缺失所需波段，缺失则退出
    if len(file_VI004) == 0:
        print('The data file of band VI004 at the input date and time is missing!')
        return False
    if len(file_VI005) == 0:
        print('The data file of band VI005 at the input date and time is missing!')
        return False
    if len(file_VI006) == 0:
        print('The data file of band VI006 at the input date and time is missing!')
        return False
    if len(file_VI008) == 0:
        print('The data file of band VI008 at the input date and time is missing!')
        return False
    if len(file_NR013) == 0:
        print('The data file of band NR013 at the input date and time is missing!')
        return False
    if len(file_NR016) == 0:
        print('The data file of band NR016 at the input date and time is missing!')
        return False
    if len(file_SW038) == 0:
        print('The data file of band SW038 at the input date and time is missing!')
        return False
    if len(file_IR112) == 0:
        print('The data file of band IR112 at the input date and time is missing!')
        return False
    if len(file_IR123) == 0:
        print('The data file of band IR123 at the input date and time is missing!')
        return False
    # 检测是否存在temp文件夹，不存在则创建
    temp_file_path = Ins_Dir + '/temp/'
    if not os.path.exists(temp_file_path):
        os.mkdir(temp_file_path)
    TOA_1000 = np.ndarray(shape=[4, 700, 700])
    TOA_1000_2000 = np.ndarray(shape=[4, 350, 350])
    TOA_2000 = np.ndarray(shape=[2, 350, 350])
    TOA_2000_1000 = np.ndarray(shape=[2, 700, 700])
    BT_2000 = np.ndarray(shape=[3, 350, 350])
    BT_2000_1000 = np.ndarray(shape=[3, 700, 700])
    lat_1000 = []
    lon_1000 = []
    lat_2000 = []
    lon_2000 = []
    for file in file_list_need:
        print(file)
        file_path = L1B_Dir + '/' + file
        data1, data2, channel, lat_1000, lon_1000, lat_2000, lon_2000 = GK2A_AMI.GetData(file_path, temp_file_path)
        if channel == 'VI004':
            TOA_1000[0] = data1
            TOA_1000_2000[0] = data2
        elif channel == 'VI005':
            TOA_1000[1] = data1
            TOA_1000_2000[1] = data2
        elif channel == 'VI008':
            TOA_1000[3] = data1
            TOA_1000_2000[3] = data2
        elif channel == 'NR013':
            TOA_2000[0] = data2
            TOA_2000_1000[0] = data1
        elif channel == 'NR016':
            TOA_2000[1] = data2
            TOA_2000_1000[1] = data1
        elif channel == 'SW038':
            BT_2000[0] = data2
            BT_2000_1000[0] = data1
        elif channel == 'IR112':
            BT_2000[1] = data2
            BT_2000_1000[1] = data1
        elif channel == 'IR123':
            BT_2000[2] = data2
            BT_2000_1000[2] = data1
    if file_VI006 != '':
        print(file_VI006)
        TOA_1000[2], TOA_1000_2000[2] = GK2A_AMI.GetData_VI006(L1B_Dir + '/' + file_VI006, temp_file_path)
    # 删除文件夹
    # shutil.rmtree(temp_file_path)
    # 临时输出
    DataTiff(TOA_1000, temp_file_path + 'TOA_1000.tif')
    so_z, so_a, sa_z, sa_a = GK2A_AMI.GetAngle(L1B_Dir + '/' + file_list_need[-1])
    print(so_z, so_a, sa_z, sa_a)
    so_Z = so_z * np.ones((350, 350))
    so_A = so_a * np.ones((350, 350))
    sa_Z = sa_z * np.ones((350, 350))
    sa_A = sa_a * np.ones((350, 350))
    if int(hour) > 22 or int(hour) < 10:
        print('云掩膜！')
        # 云掩膜
        # 2000m
        isCloud = Get_isCloud_AMI(TOA_1000_2000, TOA_2000)
        # 1000m
        isCloud_1000 = Get_isCloud_AMI(TOA_1000, TOA_2000_1000)
        print('相对方位角')
        A = Terra_MODIS.relative_angle(sa_A, sa_Z, so_A, so_Z)
        print('大气校正！')
        # 大气校正
        SRF_TOA = np.ndarray(shape=[5, 350, 350])
        for i in range(4):
            SRF_TOA[i] = TOA_1000_2000[i]
        SRF_TOA[4] = TOA_2000[1]
        SRF_2000 = algorithm_SRF_GK2A(SRF_TOA, sa_Z, so_Z, A, isCloud, lat_2000, lon_2000, Ins_Dir, name,
                                      dateTime)
        print('植被指数')
        # 植被指数
        # NDVI
        ndvi_1000 = NDVI(TOA_1000[3], TOA_1000[2])
        ndvi_2000 = NDVI(SRF_2000[3], SRF_2000[2])
        # RVI
        rvi = RVI(SRF_2000[3], SRF_2000[2])
        # NDVI
        # SAVI
        savi = SAVI(SRF_2000[3], SRF_2000[2], 0.5)
        # EVI
        evi = EVI(SRF_2000[3], SRF_2000[2], SRF_2000[0])
        # LAI
        landcover = gdal.Open('./ndata/IGBP_lai_1000.tif').ReadAsArray()
        lai = LAI(ndvi_1000, landcover)

        # 云掩膜
        for i in range(isCloud.shape[0]):
            for j in range(isCloud.shape[1]):
                if isCloud[i][j] == 1:
                    rvi[i][j] = -999
                    ndvi_2000[i][j] = -999
                    savi[i][j] = -999
                    evi[i][j] = -999
                    lai[i][j] = -999
        # 写入
        RVI_nc(rvi, lat_2000, lon_2000, Ins_Dir, name + dateTime)
        NDVI_nc(ndvi_2000, lat_2000, lon_2000, Ins_Dir, name + dateTime)
        SAVI_nc(savi, lat_2000, lon_2000, Ins_Dir, name + dateTime)
        EVI_nc(evi, lat_2000, lon_2000, Ins_Dir, name + dateTime)
        LAI_nc(lai, lat_1000, lon_1000, Ins_Dir, name + dateTime)


def GK2A_AMI_DayProducts(Ins_Dir, Day_Dir):
    name = '(GK2A-AMI)_'
    # 获取所需日期的逐时数据文件列表
    file_list_all = os.listdir(Ins_Dir)
    file_list_day = []
    date = file_list_all[0].split(name)[1][0:10]
    for file in file_list_all:
        Time = int(file.split(name)[1][11:13])
        if Time > 22 or Time < 10:
            file_list_day.append(file)
    file_list_RVI, file_list_NDVI, file_list_SAVI, file_list_EVI, file_list_LAI = Mon.GetDataList_AVHRR_AMI(
        file_list_day)
    if file_list_RVI:
        rvi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_RVI, 'RVI', 350)
        Mon.Write_nc_Day(rvi, lat, lon, Day_Dir, 'RVI', name, date)
    if file_list_NDVI:
        ndvi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_NDVI, 'NDVI', 350)
        Mon.Write_nc_Day(ndvi, lat, lon, Day_Dir, 'NDVI', name, date)
    if file_list_SAVI:
        savi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_SAVI, 'SAVI', 350)
        Mon.Write_nc_Day(savi, lat, lon, Day_Dir, 'SAVI', name, date)
    if file_list_EVI:
        evi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_EVI, 'EVI', 350)
        Mon.Write_nc_Day(evi, lat, lon, Day_Dir, 'EVI', name, date)
    if file_list_LAI:
        lai, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_LAI, 'LAI', 700)
        Mon.Write_nc_Day(lai, lat, lon, Day_Dir, 'LAI', name, date)


def Aqua_MODIS_MonProducts(Ins_Dir, Mon_NSR, Mon_SAT, Mon_DT, Mon_Dir):
    name = '(Aqua-MODIS)_'
    # 获取所需月份的逐日数据文件列表
    file_list_all = os.listdir(Ins_Dir)
    date = file_list_all[0].split(name)[1][0:7]
    year = date[0:4]
    month = date[5:7]
    file_list_RVI, file_list_NDVI, file_list_SAVI, file_list_GVI, file_list_EVI, file_list_LAI = Mon.GetDataList(
        file_list_all)
    rvi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_RVI, 'RVI', 1400)
    Mon.Write_nc(rvi, lat, lon, Mon_Dir, 'RVI', name, date)
    ndvi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_NDVI, 'NDVI', 1400)
    Mon.Write_nc(ndvi, lat, lon, Mon_Dir, 'NDVI', name, date)
    savi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_SAVI, 'SAVI', 1400)
    Mon.Write_nc(savi, lat, lon, Mon_Dir, 'SAVI', name, date)
    gvi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_GVI, 'GVI', 1400)
    Mon.Write_nc(gvi, lat, lon, Mon_Dir, 'GVI', name, date)
    evi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_EVI, 'EVI', 1400)
    Mon.Write_nc(evi, lat, lon, Mon_Dir, 'EVI', name, date)
    lai, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_LAI, 'LAI', 1400)
    Mon.Write_nc(lai, lat, lon, Mon_Dir, 'LAI', name, date)
    IGBP = './ndata/IGBP_2020.tif'
    Mon_GPP_NPP.GPP_month(year, month, ndvi, 0.005, Mon_SAT, Mon_DT, Mon_NSR, IGBP, name, Mon_Dir)


def Terra_MODIS_MonProducts(Ins_Dir, Mon_NSR, Mon_SAT, Mon_DT, Mon_Dir):
    name = '(Terra-MODIS)_'
    # 获取所需月份的逐日数据文件列表
    file_list_all = os.listdir(Ins_Dir)
    date = file_list_all[0].split(name)[1][0:7]
    year = date[0:4]
    month = date[5:7]
    file_list_RVI, file_list_NDVI, file_list_SAVI, file_list_GVI, file_list_EVI, file_list_LAI = Mon.GetDataList(
        file_list_all)
    rvi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_RVI, 'RVI', 1400)
    Mon.Write_nc(rvi, lat, lon, Mon_Dir, 'RVI', name, date)
    ndvi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_NDVI, 'NDVI', 1400)
    Mon.Write_nc(ndvi, lat, lon, Mon_Dir, 'NDVI', name, date)
    savi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_SAVI, 'SAVI', 1400)
    Mon.Write_nc(savi, lat, lon, Mon_Dir, 'SAVI', name, date)
    gvi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_GVI, 'GVI', 1400)
    Mon.Write_nc(gvi, lat, lon, Mon_Dir, 'GVI', name, date)
    evi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_EVI, 'EVI', 1400)
    Mon.Write_nc(evi, lat, lon, Mon_Dir, 'EVI', name, date)
    lai, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_LAI, 'LAI', 1400)
    Mon.Write_nc(lai, lat, lon, Mon_Dir, 'LAI', name, date)
    IGBP = './ndata/IGBP_2020.tif'
    Mon_GPP_NPP.GPP_month(year, month, ndvi, 0.005, Mon_SAT, Mon_DT, Mon_NSR, IGBP, name, Mon_Dir)


def NOAA19_AVHRR_MonProducts(Ins_Dir, Mon_NSR, Mon_SAT, Mon_DT, Mon_Dir):
    name = '(NOAA19-AVHRR)_'
    # 获取所需月份的逐日数据文件列表
    file_list_all = os.listdir(Ins_Dir)
    date = file_list_all[0].split(name)[1][0:7]
    year = date[0:4]
    month = date[5:7]
    file_list_RVI, file_list_NDVI, file_list_SAVI, file_list_EVI, file_list_LAI = Mon.GetDataList_AVHRR_AMI(
        file_list_all)
    rvi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_RVI, 'RVI', 700)
    Mon.Write_nc(rvi, lat, lon, Mon_Dir, 'RVI', name, date)
    ndvi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_NDVI, 'NDVI', 700)
    Mon.Write_nc(ndvi, lat, lon, Mon_Dir, 'NDVI', name, date)
    savi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_SAVI, 'SAVI', 700)
    Mon.Write_nc(savi, lat, lon, Mon_Dir, 'SAVI', name, date)
    evi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_EVI, 'EVI', 700)
    Mon.Write_nc(evi, lat, lon, Mon_Dir, 'EVI', name, date)
    lai, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_LAI, 'LAI', 700)
    Mon.Write_nc(lai, lat, lon, Mon_Dir, 'LAI', name, date)
    IGBP = './ndata/IGBP_2020.tif'
    Mon_GPP_NPP.GPP_month(year, month, ndvi, 0.01, Mon_SAT, Mon_DT, Mon_NSR, IGBP, name, Mon_Dir)


def GK2A_AMI_MonProducts(Ins_Dir, Mon_NSR, Mon_SAT, Mon_DT, Mon_Dir):
    name = '(GK2A-AMI)_'
    # 获取所需月份的逐日数据文件列表
    file_list_all = os.listdir(Ins_Dir)
    date = file_list_all[0].split(name)[1][0:7]
    year = date[0:4]
    month = date[5:7]
    file_list_RVI, file_list_NDVI, file_list_SAVI, file_list_EVI, file_list_LAI = Mon.GetDataList_AVHRR_AMI(
        file_list_all)
    rvi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_RVI, 'RVI', 350)
    Mon.Write_nc(rvi, lat, lon, Mon_Dir, 'RVI', name, date)
    ndvi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_NDVI, 'NDVI', 350)
    Mon.Write_nc(ndvi, lat, lon, Mon_Dir, 'NDVI', name, date)
    savi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_SAVI, 'SAVI', 350)
    Mon.Write_nc(savi, lat, lon, Mon_Dir, 'SAVI', name, date)
    evi, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_EVI, 'EVI', 350)
    Mon.Write_nc(evi, lat, lon, Mon_Dir, 'EVI', name, date)
    lai, lat, lon = Mon.CalData(Ins_Dir + '/', file_list_LAI, 'LAI', 700)
    Mon.Write_nc(lai, lat, lon, Mon_Dir, 'LAI', name, date)
    IGBP = './ndata/IGBP_2020.tif'
    Mon_GPP_NPP.GPP_month(year, month, ndvi, 0.02, Mon_SAT, Mon_DT, Mon_NSR, IGBP, name, Mon_Dir)


def Aqua_MODIS_YearProducts(Mon_Dir, Year_Dir):
    name = '(Aqua-MODIS)_'
    # 获取所需月份的逐日数据文件列表
    file_list_all = os.listdir(Mon_Dir)
    date = file_list_all[0].split(name)[1][0:4]
    file_list_RVI, file_list_NDVI, file_list_SAVI, file_list_GVI, file_list_EVI, file_list_LAI, file_list_GPP, file_list_NPP = Year.GetDataList(
        file_list_all)
    rvi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_RVI, 'RVI', 1400)
    Year.Write_nc(rvi, lat, lon, Year_Dir, 'RVI', name, date)
    ndvi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_NDVI, 'NDVI', 1400)
    Year.Write_nc(ndvi, lat, lon, Year_Dir, 'NDVI', name, date)
    savi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_SAVI, 'SAVI', 1400)
    Year.Write_nc(savi, lat, lon, Year_Dir, 'SAVI', name, date)
    gvi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_GVI, 'GVI', 1400)
    Year.Write_nc(gvi, lat, lon, Year_Dir, 'GVI', name, date)
    evi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_EVI, 'EVI', 1400)
    Year.Write_nc(evi, lat, lon, Year_Dir, 'EVI', name, date)
    lai, lat, lon = Year.CalData(Mon_Dir + '/', file_list_LAI, 'LAI', 1400)
    Year.Write_nc(lai, lat, lon, Year_Dir, 'LAI', name, date)
    GPP, lat, lon = Year.CalData(Mon_Dir + '/', file_list_GPP, 'GPP', 1400)
    Year.Write_nc(GPP, lat, lon, Year_Dir, 'GPP', name, date)
    NPP, lat, lon = Year.CalData(Mon_Dir + '/', file_list_NPP, 'NPP', 1400)
    Year.Write_nc(NPP, lat, lon, Year_Dir, 'NPP', name, date)


def Terra_MODIS_YearProducts(Mon_Dir, Year_Dir):
    name = '(Terra-MODIS)_'
    # 获取所需月份的逐日数据文件列表
    file_list_all = os.listdir(Mon_Dir)
    date = file_list_all[0].split(name)[1][0:4]
    file_list_RVI, file_list_NDVI, file_list_SAVI, file_list_GVI, file_list_EVI, file_list_LAI, file_list_GPP, file_list_NPP = Year.GetDataList(
        file_list_all)
    rvi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_RVI, 'RVI', 1400)
    Year.Write_nc(rvi, lat, lon, Year_Dir, 'RVI', name, date)
    ndvi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_NDVI, 'NDVI', 1400)
    Year.Write_nc(ndvi, lat, lon, Year_Dir, 'NDVI', name, date)
    savi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_SAVI, 'SAVI', 1400)
    Year.Write_nc(savi, lat, lon, Year_Dir, 'SAVI', name, date)
    gvi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_GVI, 'GVI', 1400)
    Year.Write_nc(gvi, lat, lon, Year_Dir, 'GVI', name, date)
    evi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_EVI, 'EVI', 1400)
    Year.Write_nc(evi, lat, lon, Year_Dir, 'EVI', name, date)
    lai, lat, lon = Year.CalData(Mon_Dir + '/', file_list_LAI, 'LAI', 1400)
    Year.Write_nc(lai, lat, lon, Year_Dir, 'LAI', name, date)
    GPP, lat, lon = Year.CalData(Mon_Dir + '/', file_list_GPP, 'GPP', 1400)
    Year.Write_nc(GPP, lat, lon, Year_Dir, 'GPP', name, date)
    NPP, lat, lon = Year.CalData(Mon_Dir + '/', file_list_NPP, 'NPP', 1400)
    Year.Write_nc(NPP, lat, lon, Year_Dir, 'NPP', name, date)


def NOAA19_AVHRR_YearProducts(Mon_Dir, Year_Dir):
    name = '(NOAA19-AVHRR)_'
    # 获取所需月份的逐日数据文件列表
    file_list_all = os.listdir(Mon_Dir)
    date = file_list_all[0].split(name)[1][0:4]
    file_list_RVI, file_list_NDVI, file_list_SAVI, file_list_EVI, file_list_LAI, file_list_GPP, file_list_NPP = Year.GetDataList_AVHRR_AMI(
        file_list_all)
    rvi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_RVI, 'RVI', 700)
    Year.Write_nc(rvi, lat, lon, Year_Dir, 'RVI', name, date)
    ndvi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_NDVI, 'NDVI', 700)
    Year.Write_nc(ndvi, lat, lon, Year_Dir, 'NDVI', name, date)
    savi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_SAVI, 'SAVI', 700)
    Year.Write_nc(savi, lat, lon, Year_Dir, 'SAVI', name, date)
    evi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_EVI, 'EVI', 700)
    Year.Write_nc(evi, lat, lon, Year_Dir, 'EVI', name, date)
    lai, lat, lon = Year.CalData(Mon_Dir + '/', file_list_LAI, 'LAI', 700)
    Year.Write_nc(lai, lat, lon, Year_Dir, 'LAI', name, date)
    GPP, lat, lon = Year.CalData(Mon_Dir + '/', file_list_GPP, 'GPP', 700)
    Year.Write_nc(GPP, lat, lon, Year_Dir, 'GPP', name, date)
    NPP, lat, lon = Year.CalData(Mon_Dir + '/', file_list_NPP, 'NPP', 700)
    Year.Write_nc(NPP, lat, lon, Year_Dir, 'NPP', name, date)


def GK2A_AMI_YearProducts(Mon_Dir, Year_Dir):
    name = '(GK2A-AMI)_'
    # 获取所需月份的逐日数据文件列表
    file_list_all = os.listdir(Mon_Dir)
    date = file_list_all[0].split(name)[1][0:4]
    file_list_RVI, file_list_NDVI, file_list_SAVI, file_list_EVI, file_list_LAI, file_list_GPP, file_list_NPP = Year.GetDataList_AVHRR_AMI(
        file_list_all)
    rvi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_RVI, 'RVI', 350)
    Year.Write_nc(rvi, lat, lon, Year_Dir, 'RVI', name, date)
    ndvi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_NDVI, 'NDVI', 350)
    Year.Write_nc(ndvi, lat, lon, Year_Dir, 'NDVI', name, date)
    savi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_SAVI, 'SAVI', 350)
    Year.Write_nc(savi, lat, lon, Year_Dir, 'SAVI', name, date)
    evi, lat, lon = Year.CalData(Mon_Dir + '/', file_list_EVI, 'EVI', 350)
    Year.Write_nc(evi, lat, lon, Year_Dir, 'EVI', name, date)
    lai, lat, lon = Year.CalData(Mon_Dir + '/', file_list_LAI, 'LAI', 700)
    Year.Write_nc(lai, lat, lon, Year_Dir, 'LAI', name, date)
    GPP, lat, lon = Year.CalData(Mon_Dir + '/', file_list_GPP, 'GPP', 350)
    Year.Write_nc(GPP, lat, lon, Year_Dir, 'GPP', name, date)
    NPP, lat, lon = Year.CalData(Mon_Dir + '/', file_list_NPP, 'NPP', 350)
    Year.Write_nc(NPP, lat, lon, Year_Dir, 'NPP', name, date)
