import numpy as np
import math
from osgeo import gdal
import os
import struct
from tqdm import tqdm
import inspect
import re
import netCDF4 as nc
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
from osgeo import gdal
from pyhdf.SD import *


# 获取日间文件列表
def GetDaytimeFiles(file_list_inDate, Date_str):
    file_list_day = list()
    for file in file_list_inDate:
        hour = int(file.split(Date_str)[1][1:3])
        if 6 < hour < 18:
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
    vrtFile = gdal.Translate(VRT_path + 'NOAA19.vrt', data, format='vrt')


# 修改VRT文件
def modify_VRT(VRT_path):
    tree = ET.parse(VRT_path + 'NOAA19.vrt')
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
    tree.write(VRT_path + 'NOAA19.vrt')


# 几何校正
def GeoCorrection(TIF, VRT, Nodata):
    correctData = gdal.Warp(TIF, VRT, dstSRS="EPSG:4490", format='GTiff',
                            xRes=0.01, yRes=0.01, outputBounds=(115, 22,
                                                                122, 29),
                            dstNodata=Nodata, srcNodata=Nodata,
                            geoloc=True, resampleAlg=gdal.GRIORA_Bilinear)
    data = correctData.ReadAsArray()
    correctData = None
    return data


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


# 单个文件经纬度范围判断
def isNeed(latitude, longitude):
    lat_min = min(np.hstack(latitude))
    lat_max = max(np.hstack(latitude))
    lon_min = min(np.hstack(longitude))
    lon_max = max(np.hstack(longitude))
    if lat_min > 29 or lat_max < 22 or lon_min > 122 or lon_max < 115:
        return False
    else:
        return True


# 单个文件数据预处理
def GetData(Ref, emissive, sa_z, so_z, A, latitude, longitude, filePath_TIF, hour):
    GetLocationTiff(longitude, latitude, filePath_TIF)
    # Ref(1,2,3b)
    GetDataTiff(Ref, filePath_TIF)
    Create_VRT(filePath_TIF + 'data.tif', filePath_TIF)
    modify_VRT(filePath_TIF)
    output_file_path = filePath_TIF + '/Output/'
    if not os.path.exists(output_file_path):
        os.mkdir(output_file_path)
    GeoCorrection(output_file_path + 'NOAA19_Ref_' + hour + '.tif', filePath_TIF + 'NOAA19.vrt', -999)
    # emissive(4,5)
    GetDataTiff(emissive, filePath_TIF)
    Create_VRT(filePath_TIF + 'data.tif', filePath_TIF)
    modify_VRT(filePath_TIF)
    output_file_path = filePath_TIF + '/Output/'
    GeoCorrection(output_file_path + 'NOAA19_emissive_' + hour + '.tif', filePath_TIF + 'NOAA19.vrt', -999)
    angle = np.ndarray(shape=[3, so_z.shape[0], so_z.shape[1]])
    angle[0] = sa_z
    angle[1] = so_z
    angle[2] = A
    # angle
    GetDataTiff(angle, filePath_TIF)
    Create_VRT(filePath_TIF + 'data.tif', filePath_TIF)
    modify_VRT(filePath_TIF)
    GeoCorrection(output_file_path + 'NOAA19_angle_' + hour + '.tif', filePath_TIF + 'NOAA19.vrt', 25500)
    TOA = gdal.Open(output_file_path + 'NOAA19_Ref_' + hour + '.tif').ReadAsArray()
    BT = gdal.Open(output_file_path + 'NOAA19_emissive_' + hour + '.tif').ReadAsArray()
    lat, lon = GetGeolocation(output_file_path + 'NOAA19_Ref_' + hour + '.tif')
    angle_data = gdal.Open(output_file_path + 'NOAA19_angle_' + hour + '.tif').ReadAsArray()
    sa_Z = angle_data[0]
    so_Z = angle_data[1]
    AA = angle_data[2]
    return TOA, BT, lat, lon, sa_Z, so_Z, AA


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


def TOA_nc(ref, emissive, sa_z, so_z, A, lat_nc, lon_nc, isCloud, out_dir, date):
    out_path = out_dir + '/' + 'TOA(NOAA19-AVHRR)_' + date + '.nc'
    f_w = nc.Dataset(out_path, 'w', format='NETCDF4')
    f_w.createDimension('lat', len(lat_nc))
    f_w.createDimension('lon', len(lon_nc))
    f_w.createDimension('size1', ref.shape[0])
    f_w.createDimension('size2', emissive.shape[0])
    f_w.createDimension('size3', 1)

    group1 = f_w.createGroup('TOA_1000')
    group1.createVariable('lat', np.float32, ('lat'))
    group1.variables['lat'][:] = lat_nc
    group1.createVariable('lon', np.float32, ('lon'))
    group1.variables['lon'][:] = lon_nc
    reflect = group1.createVariable('Top_of_Atmosphere_Reflectance_1000', np.float32, ('size1', 'lat', 'lon'))
    group1.variables['Top_of_Atmosphere_Reflectance_1000'][:] = ref
    reflect.time = date
    reflect.band_names = [1, 2, '3b']

    group2 = f_w.createGroup('BT_1000')
    group2.createVariable('lat', np.float32, ('lat'))
    group2.variables['lat'][:] = lat_nc
    group2.createVariable('lon', np.float32, ('lon'))
    group2.variables['lon'][:] = lon_nc
    emi = group2.createVariable('Brightness_Temperature_1000', np.float32, ('size2', 'lat', 'lon'))
    group2.variables['Brightness_Temperature_1000'][:] = emissive
    emi.time = date
    emi.band_names = [4, 5]

    group3 = f_w.createGroup('isCloud')
    group3.createVariable('lat', np.float32, ('lat'))
    group3.variables['lat'][:] = lat_nc
    group3.createVariable('lon', np.float32, ('lon'))
    group3.variables['lon'][:] = lon_nc
    TOA_1000 = group3.createVariable('isCloud', np.float32, ('size3', 'lat', 'lon'))
    group3.variables['isCloud'][:] = isCloud
    TOA_1000.time = date
    TOA_1000.band_names = ['isCloud']

    # 角度
    group5 = f_w.createGroup('angle')

    group5.createVariable('lat', np.float32, ('lat'))
    group5.variables['lat'][:] = lat_nc

    group5.createVariable('lon', np.float32, ('lon'))
    group5.variables['lon'][:] = lon_nc

    soz = group5.createVariable('SolarZenith', np.float32,
                                ('size3', 'lat', 'lon'))
    group5.variables['SolarZenith'][:] = so_z
    soz.time = date
    soz.band_names = ['SolarZenith']

    saz = group5.createVariable('SensorZenith', np.float32,
                                ('size3', 'lat', 'lon'))
    group5.variables['SensorZenith'][:] = sa_z
    saz.time = date
    saz.band_names = ['SensorZenith']

    a = group5.createVariable('RelativeAzimuth', np.float32,
                              ('size3', 'lat', 'lon'))
    group5.variables['RelativeAzimuth'][:] = A
    a.time = date
    a.band_names = ['RelativeAzimuth']

    f_w.close()


def expand_data(data):
    data_expand = np.zeros(2048)
    for i in range(2048):
        if i < 24:
            x1 = data[0]
            x2 = data[1]
            k = (x2 - x1)
            decimal = 0 - float((i - 24) / 40)
            data_expand[i] = k * (-decimal) + x1

        if i >= 24 and i < 2024:
            index = float((i - 24) / 40)

            decimal = index - int(index)
            data_expand[i] = data[int(index)] * (1 - decimal) + data[int(index) + 1] * decimal

        if i >= 2024:
            x1 = data[49]
            x2 = data[50]
            k = (x2 - x1)
            decimal = float((i - 24) / 40) - 50
            data_expand[i] = k * decimal + x2

    data_expand = data_expand.tolist()
    return data_expand


def create_tif(data, out_dir, name):
    data_tif = gdal.GetDriverByName("GTiff")

    # 暂时存储计算中产生的文件，包括波段数据、latitude数据和longitude数据
    data_filename = str(out_dir) + '/' + str(name) + '.tif'

    data_dataset = data_tif.Create(data_filename, data.shape[1], data.shape[0], 1, gdal.GDT_Float32)

    band0 = data_dataset.GetRasterBand(1)

    band0.WriteArray(data)

    return data_filename


def varname(p):
    for line in inspect.getframeinfo(inspect.currentframe().f_back)[3]:
        m = re.search(r'\bvarname\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)', line)
        if m:
            return m.group(1)


# 生成xml文件
def xml_file(data, lat, lon, xml, x, y):
    dom = minidom.getDOMImplementation().createDocument(None, 'VRTDataset', None)
    root = dom.documentElement
    root.setAttribute('rasterXSize', str(x))
    root.setAttribute('rasterYSize', str(y))

    element = dom.createElement('Metadata')
    element.setAttribute('domain', 'GEOLOCATION')
    root.appendChild(element)

    file = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AXIS["Lat",NORTH],AXIS["Long",EAST],AUTHORITY["EPSG","4326"]]'
    element_1 = dom.createElement('MDI')
    element_1.appendChild(dom.createTextNode(file))
    element_1.setAttribute('key', 'SRS')
    element.appendChild(element_1)

    element_2 = dom.createElement('MDI')
    element_2.appendChild(dom.createTextNode(lon))
    element_2.setAttribute('key', 'X_DATASET')
    element.appendChild(element_2)

    element_3 = dom.createElement('MDI')
    element_3.appendChild(dom.createTextNode('1'))
    element_3.setAttribute('key', 'X_BAND')
    element.appendChild(element_3)

    element_4 = dom.createElement('MDI')
    element_4.appendChild(dom.createTextNode('0'))
    element_4.setAttribute('key', 'PIXEL_OFFSET')
    element.appendChild(element_4)

    element_5 = dom.createElement('MDI')
    element_5.appendChild(dom.createTextNode('1'))
    element_5.setAttribute('key', 'PIXEL_STEP')
    element.appendChild(element_5)

    element_6 = dom.createElement('MDI')
    element_6.appendChild(dom.createTextNode(lat))
    element_6.setAttribute('key', 'Y_DATASET')
    element.appendChild(element_6)

    element_7 = dom.createElement('MDI')
    element_7.appendChild(dom.createTextNode('1'))
    element_7.setAttribute('key', 'Y_BAND')
    element.appendChild(element_7)

    element_8 = dom.createElement('MDI')
    element_8.appendChild(dom.createTextNode('0'))
    element_8.setAttribute('key', 'LINE_OFFSET')
    element.appendChild(element_8)

    element_9 = dom.createElement('MDI')
    element_9.appendChild(dom.createTextNode('1'))
    element_9.setAttribute('key', 'LINE_STEP')
    element.appendChild(element_9)

    element2 = dom.createElement('VRTRasterBand')
    element2.setAttribute('band', '1')
    element2.setAttribute('dataType', 'float32')
    root.appendChild(element2)

    element2_1 = dom.createElement('ColorInterp')
    element2_1.appendChild(dom.createTextNode('Gray'))
    element2.appendChild(element2_1)

    element2_2 = dom.createElement('NoDataValue')
    element2_2.appendChild(dom.createTextNode(''))
    element2.appendChild(element2_2)

    element2_3 = dom.createElement('SimpleSource')
    element2_3.appendChild(dom.createTextNode(''))
    element2.appendChild(element2_3)

    element2_3_1 = dom.createElement('SourceFilename')
    element2_3_1.appendChild(dom.createTextNode(data))
    element2_3_1.setAttribute('relativeToVRT', '0')
    element2_3.appendChild(element2_3_1)

    element2_3_2 = dom.createElement('SourceBand')
    element2_3_2.appendChild(dom.createTextNode('1'))
    element2_3.appendChild(element2_3_2)

    # element2_3_3 = dom.createElement('SrcRect')
    # element2_3_3.appendChild(dom.createTextNode(''))
    # element2_3_3.setAttribute('xOff', '0')
    # element2_3_3.setAttribute('yOff', '0')
    # element2_3_3.setAttribute('xSize', str(x))
    # element2_3_3.setAttribute('ySize', str(y))
    # element2_3.appendChild(element2_3_3)
    #
    # element2_3_4 = dom.createElement('DstRect')
    # element2_3_4.appendChild(dom.createTextNode(''))
    # element2_3_4.setAttribute('xOff', '0')
    # element2_3_4.setAttribute('yOff', '0')
    # element2_3_4.setAttribute('xSize', str(x))
    # element2_3_4.setAttribute('ySize', str(y))
    # element2_3.appendChild(element2_3_4)
    # 保存文件
    with open(xml, 'w+') as f:
        dom.writexml(f, addindent='\t', newl='\n')
        all_the_lines = f.readlines()
    return f


# xml文件转vrt
def renaming(file):
    """修改后缀"""

    ext = os.path.splitext(file)  # 将文件名路径与后缀名分开
    if ext[1] == '.xml':  # 文件名：ext[0]
        new_name = ext[0] + '.vrt'  # 文件后缀：ext[1]
        if os.path.exists(new_name):
            os.remove(new_name)
        os.rename(file, new_name)  # tree()已切换工作地址，直接替换后缀


def swath_georeference_vrt(out_dir, name, outputbounds, xres, yres):
    # 进行几何校正xRes=0.01, yRes=0.01

    out_filename = out_dir + '/' + 'after' + name + '.tif'
    vrt_path = out_dir + '/' + 'vrt.vrt'
    ds = gdal.Warp(out_filename, vrt_path, format='GTiff', outputBounds=outputbounds, dstSRS='EPSG:4326', xRes=xres,
                   yRes=yres, geoloc=True, resampleAlg=gdal.gdalconst.GRA_Bilinear)
    os.remove(vrt_path)

    return ds, out_filename


def read_lat_lon(ds):
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    geotransform = ds.GetGeoTransform()
    originX = geotransform[0]  # lon
    originY = geotransform[3]  # lat
    pixelWidth = geotransform[1]  # lon
    pixelHeight = geotransform[5]  # lat
    # print(originX, originY, pixelWidth, pixelHeight)
    lat = np.zeros(rows)
    lon = np.zeros(cols)

    for i in range(rows):
        lat[i] = originY + pixelHeight * i

    for i in range(cols):
        lon[i] = originX + pixelWidth * i

    return lat, lon


# 下面为读取noaa所需的一系列函数
def noaa_col(data_all):
    col1 = data_all[10:15]
    col2 = data_all[25:30]
    col3 = data_all[40:45]

    col1[0] = col1[0] / 10 ** 10
    col2[0] = col2[0] / 10 ** 10
    col3[0] = col3[0] / 10 ** 10

    col1[1] = col1[1] / 10 ** 7
    col2[1] = col2[1] / 10 ** 7
    col3[1] = col3[1] / 10 ** 7

    col1[2] = col1[2] / 10 ** 10
    col2[2] = col2[2] / 10 ** 10
    col3[2] = col3[2] / 10 ** 10

    col1[3] = col1[3] / 10 ** 7
    col2[3] = col2[3] / 10 ** 7
    col3[3] = col3[3] / 10 ** 7

    col1[4] = col1[4]
    col2[4] = col2[4]
    col3[4] = col3[4]

    col = list()
    col.append(col1)
    col.append(col2)
    col.append(col3)

    col = np.array(col).astype(float)

    return col


def noaa_ir(ir_all):
    ir1 = ir_all[0:3]
    ir2 = ir_all[6:9]
    ir3 = ir_all[12:15]

    ir1[0] = ir1[0] / 10 ** 9
    ir1[1] = ir1[1] / 10 ** 6
    ir1[2] = ir1[2] / 10 ** 6

    ir2[0] = ir2[0] / 10 ** 9
    ir2[1] = ir2[1] / 10 ** 6
    ir2[2] = ir2[2] / 10 ** 6

    ir3[0] = ir3[0] / 10 ** 9
    ir3[1] = ir3[1] / 10 ** 6
    ir3[2] = ir3[2] / 10 ** 6

    ir = list()

    ir.append(ir1)
    ir.append(ir2)
    ir.append(ir3)

    ir = np.array(ir).astype(float)

    return ir


def noaa_angle(angle_all):
    so = angle_all[0::3]
    sa = angle_all[1::3]
    a = angle_all[2::3]

    so = np.array(so).astype(float) / 100
    sa = np.array(sa).astype(float) / 100
    a = np.array(a).astype(float) / 100
    return so, sa, a


def geo_avhrr(geo):
    lat = geo[0::2]
    lon = geo[1::2]

    lat = np.array(lat).astype(float) / 10 ** 4
    lon = np.array(lon).astype(float) / 10 ** 4

    return lat, lon


def band_avhrr(band):
    band1 = band[0::5]
    band2 = band[1::5]
    band3 = band[2::5]
    band4 = band[3::5]
    band5 = band[4::5]

    band1 = np.array(band1).astype(float)
    band2 = np.array(band2).astype(float)
    band3 = np.array(band3).astype(float)
    band4 = np.array(band4).astype(float)
    band5 = np.array(band5).astype(float)

    return band1, band2, band3, band4, band5


def noaa_read(file_path):
    ds = open(file_path, "rb")
    size = os.path.getsize(file_path)
    skip = 22016
    so_all = list()
    sa_all = list()
    a_all = list()
    lat_all = list()
    lon_all = list()
    band1_all = list()
    band2_all = list()
    band3_all = list()
    band4_all = list()
    band5_all = list()
    print('开始逐行读取数据！')
    for num in tqdm(range(int(size / 22016))):
        if num == 0:
            # year
            ds.seek(84, 0)
            data = ds.read(2)
            r = float(struct.unpack('H', data[0:2])[0])
            print(r)
            # month
            ds.seek(86, 0)
            data = ds.read(2)
            r = float(struct.unpack('H', data[0:2])[0])
            print(r)
            # day
            ds.seek(88, 0)
            data = ds.read(4)
            r = float(struct.unpack('I', data[0:4])[0])
            print(r / 3600000)

            ds.seek(292, 0)
            data = ds.read(12)
            ch4 = list()
            for i in range(3):
                r = float(struct.unpack('i', data[i * 4:i * 4 + 4])[0])
                ch4.append(r)

            ds.seek(304, 0)
            data = ds.read(12)
            ch5 = list()
            for i in range(3):
                r = float(struct.unpack('i', data[i * 4:i * 4 + 4])[0])
                ch5.append(r)

        else:
            ds.seek(num * skip + 48, 0)
            data = ds.read(180)
            col_all = list()
            for i in range(45):
                col_s = float(struct.unpack('i', data[4 * i:4 * i + 4])[0])
                col_all.append(col_s)
            col = noaa_col(col_all)

            # 读取红外波段定标系数
            ds.seek(num * skip + 228, 0)
            data = ds.read(72)
            ir_all = list()
            for i in range(18):
                ir_s = (struct.unpack('i', data[4 * i:4 * i + 4])[0])
                ir_all.append(ir_s)

            ir = noaa_ir(ir_all)

            # 读取天顶角、方位角
            ds.seek(num * skip + 328, 0)
            data = ds.read(306)
            angle_all = list()
            for i in range(153):
                angle_s = float(struct.unpack('h', data[2 * i:2 * i + 2])[0])
                angle_all.append(angle_s)

            so, sa, a = noaa_angle(angle_all)
            so_expand = expand_data(so)
            sa_expand = expand_data(sa)
            a_expand = expand_data(a)

            so_all.append(so_expand)
            sa_all.append(sa_expand)
            a_all.append(a_expand)

            # 读取经纬度
            ds.seek(num * skip + 640, 0)
            data = ds.read(408)
            geo_all = list()
            for i in range(102):
                geo_s = float(struct.unpack('i', data[4 * i:4 * i + 4])[0])
                geo_all.append(geo_s)

            lat, lon = geo_avhrr(geo_all)

            lat_expand = expand_data(lat)
            lon_expand = expand_data(lon)

            lat_all.append(lat_expand)
            lon_all.append(lon_expand)

            # 读取波段数据
            ds.seek(num * skip + 1264, 0)
            data = ds.read(20480)
            band_all = list()
            for i in range(10240):
                band_s = float(struct.unpack('h', data[2 * i:2 * i + 2])[0])
                band_all.append(band_s)

            band1, band2, band3, band4, band5 = band_avhrr(band_all)

            c1 = 1.1910427 / 10 ** 5
            c2 = 1.4387752
            # print(band1[0:10])
            # print(band2[0:10])
            # print(band3[0:10])
            # print(band5[0:10])
            # print(band3[0:10])
            # print(col)
            # 进行辐射校正
            for i in range(band1.shape[0]):
                if band1[i] > col[0][4]:
                    band1[i] = band1[i] * col[0][2] + col[0][3]
                else:
                    band1[i] = band1[i] * col[0][0] + col[0][1]

                if band2[i] > col[1][4]:
                    band2[i] = band2[i] * col[1][2] + col[1][3]
                else:
                    band2[i] = band2[i] * col[1][0] + col[1][1]

                if band3[i] > col[2][4]:
                    band3[i] = band3[i] * col[2][2] + col[2][3]
                else:
                    band3[i] = band3[i] * col[2][0] + col[2][1]

                band4[i] = band4[i] ** 2 * ir[1][0] + band4[i] * ir[1][1] + ir[1][2]
                wave = ch4[0] / 1000
                band4[i] = c2 * wave / (math.log(math.fabs(1 + (c1 * wave * wave * wave / band4[i]))))
                a = ch4[1] / 10 ** 5
                b = ch4[2] / 10 ** 6
                band4[i] = (band4[i] - a) / b

                band5[i] = band5[i] ** 2 * ir[2][0] + band5[i] * ir[2][1] + ir[2][2]
                wave = ch5[0] / 1000
                band5[i] = c2 * wave / (math.log(math.fabs(1 + (c1 * wave * wave * wave / band5[i]))))
                a = ch5[1] / 10 ** 5
                b = ch5[2] / 10 ** 6
                band5[i] = (band5[i] - a) / b

            band1_all.append(band1)
            band2_all.append(band2)
            band3_all.append(band3)
            band4_all.append(band4)
            band5_all.append(band5)

    band1_all = np.array(band1_all).astype(float)
    band2_all = np.array(band2_all).astype(float)
    band3_all = np.array(band3_all).astype(float)
    band4_all = np.array(band4_all).astype(float)
    band5_all = np.array(band5_all).astype(float)

    ref = list()
    ref.append(band1_all / 100)
    ref.append(band2_all / 100)
    ref.append(band3_all / 100)

    irt = list()
    irt.append(band4_all)
    irt.append(band5_all)

    ref = np.array(ref).astype(float)
    irt = np.array(irt).astype(float)

    sa_all = np.array(sa_all).astype(float)
    so_all = np.array(so_all).astype(float)
    a_all = np.array(a_all).astype(float)
    lat_all = np.array(lat_all).astype(float)
    lon_all = np.array(lon_all).astype(float)
    print(ref[0])
    for i in range(ref.shape[1]):
        for j in range(ref.shape[2]):
            ref[0][i][j] = ref[0][i][j] / math.cos(so_all[i][j] * math.pi / 180)
            ref[1][i][j] = ref[1][i][j] / math.cos(so_all[i][j] * math.pi / 180)
    return ref, irt, sa_all, so_all, a_all, lat_all, lon_all
