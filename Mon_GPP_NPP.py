import datetime
from netCDF4 import Dataset
import numpy as np
import sys
from osgeo import gdal
import netCDF4 as nc
import numpy as np
import sys
from osgeo import gdal
from PIL import Image
from scipy import ndimage


def extract_pixel_at_PxPy_from_tif(src_file, px, py):
    # extract point from geotif file
    # src_file：tif文件名、px，py：坐标值
    src_ds = gdal.Open(src_file)
    gt = src_ds.GetGeoTransform()
    rb = src_ds.GetRasterBand(1)
    intval = rb.ReadAsArray(px, py, 1, 1)
    return (intval[0][0])


def extract_pixel_at_PxPy_from_netcdf(src_file, px, py, var_name):
    # extract point from netcdf file
    # src_file：nc文件名、px，py：坐标值,var_name:变量名
    Data_nc = Dataset(src_file)
    Data = Data_nc[var_name]
    lats = Data_nc.variables['latitude'][:]
    lons = Data_nc.variables['longitude'][:]

    indx = 0
    for i in range(len(lats)):
        if (lats[i] < py):
            break
    indx = i

    indy = 0
    for j in range(len(lons)):
        if (lons[j] > px):
            break
    indy = j

    # print(indx,indy)
    return (Data[0, indx, indy])


def y_m(year, month):
    # 返回该年月的天数
    # year：该日期的年份, month：该日期的月份
    if month > 12 or month <= 0:
        print("月份有误！")
        return -1
    if month == 2:
        if year % 4 == 0 and year % 100 != 0 or year % 400 == 0:
            return 29
        else:
            return 28
    elif month in (4, 6, 9, 11):
        return 30
    else:
        return 31


def GPP_NPP_pixel(year, month, T2, DT, SSR, ndvi, IGBP_type):
    # 返回某日期当月的某个经纬度的GPP和NPP逐月值
    # year：该日期的年份, month：该日期的月份
    # T2：2m气温像素值, DT：露点温度像素值, SSR：净太阳辐射像素值,IGBP_type:像素地类代码
    # T2 in K; DT in K, SSR in J m-2
    # 地类字典
    a_b_dic = {2: [1.1505, -0.1294], 4: [1.1505, -0.1294], 5: [1.1505, -0.1294],
               8: [1.0705, -0.0141], 9: [1.0705, -0.0141], 10: [1.1225, -0.0221],
               12: [1.0344, -0.0268], 14: [1.1225, -0.0221]
               }  # 根据ndvi计算fpar的斜率a和截距b

    igbp = list([2, 4, 5, 8, 9, 10, 12, 14])
    epsilon_max = [1.268, 1.165, 1.051, 1.239, 1.206, 0.86, 1.044, 1.044]
    Tminmin = [-8, -6, -8, -8, -8, -8, -8, -8]
    Tminmax = [9.09, 9.94, 11.39, 11.39, 11.39, 10.44, 12.02, 12.02]
    VPDmin = [800, 650, 650, 650, 650, 650, 650, 650]
    VPDmax = [3100, 1650, 3200, 3200, 3100, 2300, 4300, 4300]
    if (ndvi < 0):
        GPP_1pixel = -999
        NPP_1pixel = -999
    elif (IGBP_type in igbp):
        t2m = T2 - 273.15  # 转成℃的单位
        d2m = DT - 273.15  # 转成℃的单位
        # 计算得出VPD
        SVP = 6.112 * np.exp((17.67 * t2m) / (t2m + 243.5)) * 100  # 单位为Pa Tetens 公式
        Rh = 100.0 * (np.exp((17.625 * d2m) / (243.03 + d2m)) / np.exp((17.625 * t2m) / (243.04 + t2m)))
        VPD = SVP * (1 - Rh / 100.0)
        # print(SVP, Rh, VPD)

        SR = SSR * 1e-6  # J m-2 to MJ m-2

        indx = igbp.index(IGBP_type)

        # calculating Tscalar
        Tscalar = (t2m - Tminmin[indx]) / (Tminmax[indx] - Tminmin[indx])
        if (Tscalar > 1):
            Tscalar = 1
        if (Tscalar < 0):
            Tscalar = 0

        # calculating VPDscalar
        VPDscalar = (VPDmax[indx] - VPD) / (VPDmax[indx] - VPDmin[indx])
        if (VPDscalar > 1):
            VPDscalar = 1
        if (VPDscalar < 0):
            VPDscalar = 0

        # fpar
        # 取出igbp对应的a和b
        a = a_b_dic[IGBP_type][0]
        b = a_b_dic[IGBP_type][1]
        fpar = ndvi * a + b
        if (fpar > 0.95):
            fpar = 0.95
        if (fpar < 0):
            fpar = 0

        # GPP!!!
        days = y_m(year, month)
        # unit: gC/m-2/month
        GPP_1pixel = epsilon_max[indx] * fpar * SR * Tscalar * VPDscalar * days

        # NPP
        Reco_ratio_1pixel = (t2m - 0) * (t2m - 35) / ((t2m - 0) * (t2m - 35) - (t2m - 20) * (t2m - 20))
        if (Reco_ratio_1pixel < 0.4):
            Reco_ratio_1pixel = 0.4
        if (Reco_ratio_1pixel > 0.7):
            Reco_ratio_1pixel = 0.7

        NPP_1pixel = GPP_1pixel * (1 - Reco_ratio_1pixel)

        # print(GPP_1pixel, NPP_1pixel, epsilon_max[indx], SR, ndvi, fpar, t2m, Tscalar, VPD, VPDscalar)

    else:
        # 空值
        GPP_1pixel = -999
        NPP_1pixel = -999

    return (GPP_1pixel, NPP_1pixel)


def creat_nc_file(filename, varname, data, lat, lon):
    # 创建nc文件
    ncfile = Dataset(filename, 'w', format='NETCDF4')
    ncfile.createDimension('lat', len(lat))  # latitude axis
    ncfile.createDimension('lon', len(lon))  # longitude axis

    ncfile.title = varname

    group = ncfile.createGroup(varname)

    latitude = group.createVariable('lat', np.float32, ('lat',))
    group.variables['lat'][:] = lat
    latitude.units = 'degrees_north'
    latitude.long_name = 'latitude'

    longitude = group.createVariable('lon', np.float32, ('lon',))
    group.variables['lon'][:] = lon
    longitude.units = 'degrees_east'
    longitude.long_name = 'longitude'

    var = group.createVariable(varname, np.float32, ('lat', 'lon'))
    group.variables[varname][:] = data
    var.units = 'gC/m2/m'  # degrees Kelvin
    var.standard_name = varname  # this is a CF standard name

    ncfile.close()


def read_ECWMF(src_file, var_name, Size):
    Data_nc = nc.Dataset(src_file)
    Data = Data_nc[var_name][0][:]
    # scale = Data_nc[var_name].scale_factor
    # offset = Data_nc[var_name].add_offset
    Data = np.array(Data)
    # Data = Data * scale + offset
    new_pil_image = ndimage.zoom(Data, (Size / Data.shape[0], Size / Data.shape[1]), order=1)
    new_pil_image = np.array(new_pil_image)

    return new_pil_image


def GPP_month(year, month, NDVIdata, ndvisize, TFilename, DTFilename, SRFilename, IGBP_file, name, outPath):
    '''计算某年某月福建省范围的GPP和NPP
    :param SRFilename:净光合辐射量文件
    :param NDVIdata:ndvi的数据数组
    :param ndvisize:ndvi的数据的栅格分辨率 度为单位
    :param TFilename:气温文件
    :param DTFilename:露点温度文件
    :param SRFilename:露点温度文件
    :param IGBP_file:GPP与NPP计算所用到的IGBP文件
    :param outPath:GPP与NPP的文件输出路径
   '''

    adfGeoTransform = [115.0, 0.05, 0.0, 29.0, 0.0, -0.05]  # 默认地理变换参数
    # 左上角地理坐标115 29
    # print(adfGeoTransform[0])
    # print(adfGeoTransform[3])

    year1 = int(year)
    month1 = int(month)

    nXSize = int(7.0 / ndvisize)  # 列数
    nYSize = int(7.0 / ndvisize)  # 行数
    progress = nXSize * nYSize
    adfGeoTransform[1] = ndvisize
    adfGeoTransform[5] = -ndvisize

    GPP = []
    NPP = []
    T2m = read_ECWMF(TFilename, 't2m', nXSize)
    DT = read_ECWMF(DTFilename, 'd2m', nXSize)
    SR = read_ECWMF(SRFilename, 'ssrc', nXSize)
    for i in range(int(nYSize)):
        for j in range(int(nXSize)):
            px = adfGeoTransform[0] + i * adfGeoTransform[1] + j * adfGeoTransform[2]
            py = adfGeoTransform[3] + i * adfGeoTransform[4] + j * adfGeoTransform[5]

            ndvi_1pixel = NDVIdata[i, j]
            vegtype_1pixel = extract_pixel_at_PxPy_from_tif(IGBP_file, px, py)  # slice nc file
            # T2m_1pixel = extract_pixel_at_PxPy_from_netcdf(TFilename, px, py, "t2m")  # slice nc file
            # DT_1pixel = extract_pixel_at_PxPy_from_netcdf(DTFilename, px, py, "d2m")  # slice nc file
            # SR_1pixel = extract_pixel_at_PxPy_from_netcdf(SRFilename, px, py, "ssrc")  # slice nc file
            T2m_1pixel = T2m[i, j]
            DT_1pixel = DT[i, j]
            SR_1pixel = SR[i, j]
            # unit: gC/m-2/months
            GPP_1p, NPP_1p = GPP_NPP_pixel(year1, month1, T2m_1pixel, DT_1pixel, SR_1pixel, ndvi_1pixel, vegtype_1pixel)
            # unit: KgC/m-2/month ,scale_factor =  0.0001
            GPP.append(GPP_1p * 0.45)
            NPP.append(NPP_1p * 0.45)
            print("\r", end="")
            p = 100 * len(GPP) / progress
            print("进度: {:.2f}%: ".format(p), end="")
            sys.stdout.flush()
    date = year + '_' + month
    # output GPP as netcdf
    GPP = np.reshape(GPP, [nYSize, nXSize])

    out_GPP_file = outPath + '/Mon-GPP' + name + date + '.nc'
    lats = list()

    for i in range(nYSize):
        py = adfGeoTransform[3] + i * adfGeoTransform[5]
        lats.append(py)

    lons = list()
    for j in range(nXSize):
        px = adfGeoTransform[0] + j * adfGeoTransform[1]
        lons.append(px)
    creat_nc_file(out_GPP_file, "GPP", GPP, lats, lons)

    # output NPP as netcdf
    NPP = np.reshape(NPP, [nYSize, nXSize])
    out_NPP_file = outPath + '/Mon-NPP' + name + date + '.nc'
    creat_nc_file(out_NPP_file, "NPP", NPP, lats, lons)
