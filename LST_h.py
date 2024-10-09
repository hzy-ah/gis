import numpy as np
import netCDF4 as nc
from osgeo import gdal
import create_data as cd

# 1.读入亮温数据T31,T32
def get_Brightness_Temperature(T31, T32):
    Ti = [T31, T32]
    return Ti


# 2.计算劈窗算法系数A0,A1,A2
# 计算植被覆盖度Pv
def get_Vegetation_Coverage(NDVI):
    NDVIv = 0.9  # 茂密植被覆盖像元的NDVI
    NDVIs = 0.15  # 完全裸土覆盖像元的NDVI
    if NDVI > NDVIv:
        Pv = 1
    elif NDVI < NDVIs:
        Pv = 0
    else:
        Pv = (NDVI - NDVIs) / (NDVIv - NDVIs)
    return Pv


# 根据像元地表覆盖类型得到第i通道发射率Evi
def get_Evi(Surface_Type):
    Ev = []
    # Ev = [0.98672, 0.98990]
    if Surface_Type == 1:
        Ev = [0.992752667, 0.992752667]
    elif Surface_Type == 2:
        Ev = [0.992752667, 0.992752667]
    elif Surface_Type == 3:
        Ev = [0.981234923, 0.981148692]
    elif Surface_Type == 4:
        Ev = [0.981234923, 0.981148692]
    elif Surface_Type == 5:
        Ev = [0.981234923, 0.981148692]
    elif Surface_Type == 6:
        Ev = [0.982393, 0.986172]
    elif Surface_Type == 7:
        Ev = [0.982393, 0.986172]
    elif Surface_Type == 8:
        Ev = [0.972976, 0.972952]
    elif Surface_Type == 9:
        Ev = [0.972976, 0.972952]
    elif Surface_Type == 10:
        Ev = [0.972976, 0.972952]
    elif Surface_Type == 11:
        Ev = [0.972976, 0.972952]
    elif Surface_Type == 12:
        Ev = [0.982393, 0.983454]
    elif Surface_Type == 13:
        Ev = [0.968944891, 0.970830127]
    elif Surface_Type == 14:
        Ev = [0.982393, 0.983454]
    elif Surface_Type == 15 or Surface_Type == 16:
        Ev = [0.992905667, 0.987093167]
    elif Surface_Type == 17:
        Ev = [0, 0]
    return Ev


# 估算热辐射相互作用校正项de
def get_de(Pv):
    de = 0
    if Pv == 0 or Pv == 1:
        de = 0
    elif 0 < Pv < 0.5:
        de = 0.0003796 * Pv
    elif 0.5 < Pv < 1:
        de = 0.0003796 * (1 - Pv)
    elif Pv == 0.5:
        de = 0.001898
    return de


# 2.1 计算地表比辐射率E31,E32
# Ei = Pv * Rv * Eiv + (1 - Pv) * Rs * Eis + de
def get_Surface_Emissivity(NDVI, Surface_Type):
    # 1.计算植被覆盖度Pv
    Pv = get_Vegetation_Coverage(NDVI)
    # 2.计算植被和裸土的辐射比率
    Rv = 0.92762 + 0.07033 * Pv
    Rs = 0.99782 + 0.08362 * Pv
    # 3.根据像元地表覆盖类型得到第i通道发射率
    Evi = get_Evi(Surface_Type)
    # 4.定义裸土在第i通道的地表比辐射率
    Esi = [0.96767, 0.97790]
    # 5.估算热辐射相互作用校正项
    de = get_de(Pv)
    # 6.计算第i波段的地表比辐射率
    Ei = [Pv * Rv * Evi[0] + (1 - Pv) * Rs * Esi[0] + de, Pv * Rv * Evi[1] + (1 - Pv) * Rs * Esi[1] + de]
    return Ei


# 2.2 计算大气透过率（夏季w=2,冬季w=1,3,5）
def get_Atmospheric_Transmittance(w):
    At = []
    At.append(1.07268 - 0.12571 * w)
    At.append(0.93821 - 0.12613 * w)
    return At


# 2.3 计算Qin劈窗算法系数A0、A1、A2
def get_Coefficient(Ei, At):
    C = []
    D = []
    A = []
    for i in range(2):
        C.append(Ei[i] * At[i])
        D.append((1 - At[i]) * (1 + (1 - Ei[i]) * At[i]))
    a31 = -64.60363
    b31 = 0.440817
    c31 = C[0]
    d31 = D[0]
    a32 = -68.72575
    b32 = 0.473453
    c32 = C[1]
    d32 = D[1]
    if (d32 * c31 - d31 * c32) == 0:
        A.append(-999)
        A.append(-999)
        A.append(-999)
    else:
        A.append((a31 * d32 * (1 - c31 - d31)) / (d32 * c31 - d31 * c32) - (a32 * d31 * (1 - c32 - d32)) / (
                d32 * c31 - d31 * c32))
        A.append(1 + d31 / (d32 * c31 - d31 * c32) + (b31 * d32 * (1 - c31 - d31)) / (d32 * c31 - d31 * c32))
        A.append(d31 / (d32 * c31 - d31 * c32) + (b32 * d31 * (1 - c31 - d31)) / (d32 * c31 - d31 * c32))

    return A


# 3.代入系数计算地表温度Ts
def get_LST(Ti, A):
    Ts = A[0] + A[1] * Ti[0] - A[2] * Ti[1]
    if A[0] == -999:
        Ts = -999
    return Ts


# 4.写入nc文件
def LST_nc(Ins_Dir, Ts, lat_1000, lon_1000, text):
    Output_path = Ins_Dir + '/' + text[0] + '_' + text[1] + '.nc'
    nc_file = nc.Dataset(Output_path, 'w', format='NETCDF4')
    nc_file.createDimension('lat', len(lat_1000))
    nc_file.createDimension('lon', len(lon_1000))

    group = nc_file.createGroup('LST')

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat_1000

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon_1000

    BT_1000 = group.createVariable('LST', np.float32,
                                   ('lat', 'lon'))
    group.variables['LST'][:] = Ts
    BT_1000.time = text[1]


def algorithm_LST(BT, NDVI, isCloud, Ins_Dir, lat_1000, lon_1000, text):
    Ts = -999 * np.ones((BT.shape[1], BT.shape[2]))
    IGBP = gdal.Open('./ndata/IGBP_LST.tif').ReadAsArray()
    for i in range(BT.shape[1]):
        for j in range(BT.shape[2]):

            if isCloud[i][j] == 1 or BT[0][i][j] == -999:
                continue

            Ti = get_Brightness_Temperature(BT[0][i][j], BT[1][i][j])

            Ei = get_Surface_Emissivity(NDVI[i][j], IGBP[i][j])
            if Ei == 0:
                Ts[i][j] = -999
                continue
            At = get_Atmospheric_Transmittance(2)
            A = get_Coefficient(Ei, At)
            Ts[i][j] = get_LST(Ti, A)
            if Ts[i][j] > 350:
                Ts[i][j] = -999
    cd.create_tif(Ts, Ins_Dir, 'lst')
    LST_nc(Ins_Dir, Ts, lat_1000, lon_1000, text)