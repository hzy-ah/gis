import netCDF4 as nc
import numpy as np
import sys
from osgeo import gdal
from PIL import Image
from scipy import ndimage

def extract_pixel_at_PxPy_from_tif(src_file, px, py):
    src_ds = gdal.Open(src_file)
    rb = src_ds.GetRasterBand(1)
    intval = rb.ReadAsArray(px, py, 1, 1)
    return (intval[0][0])


def extract_pixel_at_PxPy_from_netcdf(src_file, px, py, var_name):
    Data_nc = nc.Dataset(src_file)
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





    # print(indx,indy)(Data[0, indx, indy])
    return Data[0, indx, indy]


def y_m(year, month):
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

        GPP_1pixel = epsilon_max[indx] * Tscalar * VPDscalar * fpar * SR * days


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


def creat_nc_file(filename, varname, data, lats, lons):
    ncfile = nc.Dataset(filename, mode='w', format='NETCDF4')

    lat_dim = ncfile.createDimension('lat', len(lats))
    lon_dim = ncfile.createDimension('lon', len(lons))

    group = ncfile.createGroup(varname)
    var = group.createVariable(varname, np.float64, ('lat', 'lon'))  # note: unlimited dimension is leftmost
    var.units = 'gC/m2/m'
    var.standard_name = varname

    lat = group.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'
    lon = group.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    lon.long_name = 'longitude'

    lat[:] = lats
    lon[:] = lons
    var[:, :] = data
    ncfile.close()

def read_ECWMF(src_file, var_name, Size):
    Data_nc = nc.Dataset(src_file)
    Data = Data_nc[var_name][0][:]
    # scale = Data_nc[var_name].scale_factor
    # offset = Data_nc[var_name].add_offset
    Data = np.array(Data)
    # Data = Data * scale + offset
    new_pil_image = ndimage.zoom(Data, (Size/Data.shape[0], Size/Data.shape[1]), order=1)
    new_pil_image = np.array(new_pil_image)

    return new_pil_image

def read_IGBP(IGBP_file, Size):
    src_ds = gdal.Open(IGBP_file)
    rb = src_ds.GetRasterBand(1)
    intval = rb.ReadAsArray()
    intval = np.array(intval)

    new_pil_image = ndimage.zoom(intval, (Size / intval.shape[0], Size / intval.shape[1]), order=1)
    new_pil_image = np.array(new_pil_image)



    return new_pil_image




def GPP_month(year, month, NDVIdata, ndvisize, TFilename, DTFilename, SRFilename, IGBP_file, outPath, header):
    adfGeoTransform = [115.0, 0.05, 0.0, 29.0, 0.0, -0.05]

    year1 = int(year)
    month1 = int(month)

    nXSize = int(7.0 / ndvisize)
    nYSize = int(7.0 / ndvisize)
    progress = nXSize * nYSize
    print(nXSize)

    # 替换为真实的ndvi大小
    adfGeoTransform[1] = ndvisize
    adfGeoTransform[5] = -ndvisize

    GPP = []
    NPP = []
    T2m = read_ECWMF(TFilename, 't2m', nXSize)
    DT = read_ECWMF(DTFilename, 'd2m', nXSize)
    SR = read_ECWMF(SRFilename, 'ssrc', nXSize)
    # vegtype = read_IGBP(IGBP_file, nXSize)



    for i in range(int(nYSize)):
        for j in range(int(nXSize)):

            px = adfGeoTransform[0] + i * adfGeoTransform[1] + j * adfGeoTransform[2]
            py = adfGeoTransform[3] + i * adfGeoTransform[4] + j * adfGeoTransform[5]

            ndvi_1pixel = NDVIdata[i, j]

            vegtype_1pixel = extract_pixel_at_PxPy_from_tif(IGBP_file, px, py)
            # T2m_1pixel = extract_pixel_at_PxPy_from_netcdf(TFilename, px, py, "t2m")
            # DT_1pixel = extract_pixel_at_PxPy_from_netcdf(DTFilename, px, py, "d2m")
            # SR_1pixel = extract_pixel_at_PxPy_from_netcdf(SRFilename, px, py, "ssrc")
            # vegtype_1pixel = vegtype[i, j]
            T2m_1pixel = T2m[i, j]
            DT_1pixel = DT[i, j]
            SR_1pixel = SR[i, j]

            GPP_1p, NPP_1p = GPP_NPP_pixel(year1, month1, T2m_1pixel, DT_1pixel, SR_1pixel, ndvi_1pixel, vegtype_1pixel)

            GPP.append(GPP_1p * 0.45)
            NPP.append(NPP_1p * 0.45)
            print("\r", end="")
            p = 100 * len(GPP) / progress
            print("进度: {:.2f}%: ".format(p), end="")
            sys.stdout.flush()

    GPP = np.reshape(GPP, [nYSize, nXSize])

    out_GPP_file = outPath + '/' + header[0] + '_' + year + '_' + month + '.nc'
    lats = list()

    for i in range(nYSize):
        py = adfGeoTransform[3] + i * adfGeoTransform[5]
        lats.append(py)

    lons = list()
    for j in range(nXSize):
        px = adfGeoTransform[0] + j * adfGeoTransform[1]
        lons.append(px)

    creat_nc_file(out_GPP_file, "GPP", GPP, lats, lons)

    NPP = np.reshape(NPP, [nYSize, nXSize])
    out_NPP_file = outPath + '/' + header[1] + '_' + year + '_' + month + '.nc'
    creat_nc_file(out_NPP_file, "NPP", NPP, lats, lons)
