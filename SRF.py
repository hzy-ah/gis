import numpy as np
import netCDF4 as nc
from osgeo import gdal


# 读取查找表数据
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


def h2aod_id(h):
    if h < 1:
        aod_id = 1
    if 1 <= h <= 2:
        aod_id = 1
    if h > 2:
        aod_id = 0
    return aod_id


# 大气校正函数
# 带改进的文件，读取查找表nc文件的元数据，直接获取
def atmos_correct(data, my_Altitude, sa_z, so_z, A, Pa, Ra, TT, cloud):
    Pa_PList = np.zeros((sa_z.shape[0], sa_z.shape[1]))
    Ra_PList = np.zeros((sa_z.shape[0], sa_z.shape[1]))
    TT_PList = np.zeros((sa_z.shape[0], sa_z.shape[1]))
    Rs = np.zeros((sa_z.shape[0], sa_z.shape[1]))

    for i in range(sa_z.shape[0]):
        for j in range(sa_z.shape[1]):
            if cloud[i][j] == 0:
                H = my_Altitude[i][j]
                my_Mu0 = so_z[i][j]
                my_Mu = sa_z[i][j]
                my_RAA = A[i][j]

                if H > 3:
                    H = 3

                if H < 0:
                    H = 0

                H_id = round(H)
                my_AOD_id = h2aod_id(H)
                my_Mu0 = 2 * round(my_Mu0 / 2)
                my_Mu = 2 * round(my_Mu / 2)
                my_RAA = 5 * round(my_RAA / 5)

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
                if Rs[i][j] < 0:
                    print(data[i][j],my_Mu,my_RAA,my_Mu0)

            if cloud[i][j] == 1:
                Rs[i][j] = -999

    return Rs


def atmos(ref, sa_z, so_z, A, table_file, cloud):
    path = './ndata/dem.tif'
    f = gdal.Open(path)
    h = f.ReadAsArray()
    h = np.array(h).astype(float)

    if cloud.shape[0] == 700:
        h = h[::2, ::2]
    elif cloud.shape[1] == 350:
        h = h[::4, ::4]

    Pa, Ra, TT, AOD, Altitude, Mu, Mu0, RAA = table_data(table_file)

    Rs = atmos_correct(ref, h, sa_z, so_z, A, Pa, Ra, TT, cloud)

    return Rs


def SRF_MODIS_nc(SRF, lat, lon, Ins_Dir, name, date):
    Output_path = Ins_Dir + '/' + 'SRF' + name + date + '.nc'
    nc_file = nc.Dataset(Output_path, 'w', format='NETCDF4')
    nc_file.createDimension('lat', len(lat))
    nc_file.createDimension('lon', len(lon))
    nc_file.createDimension('size', SRF.shape[0])

    group = nc_file.createGroup('SRF')

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon

    Surface_Reflectance = group.createVariable('Surface_Reflectance', np.float32,
                                               ('size', 'lat', 'lon'))
    group.variables['Surface_Reflectance'][:] = SRF
    Surface_Reflectance.time = date
    Surface_Reflectance.band_names = [1, 2, 3, 4, 5, 6, 7]


def SRF_NOAA_nc(SRF, lat, lon, Ins_Dir, name, date):
    Output_path = Ins_Dir + '/' + 'SRF' + name + date + '.nc'
    nc_file = nc.Dataset(Output_path, 'w', format='NETCDF4')
    nc_file.createDimension('lat', len(lat))
    nc_file.createDimension('lon', len(lon))
    nc_file.createDimension('size', SRF.shape[0])

    group = nc_file.createGroup('SRF')

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon

    Surface_Reflectance = group.createVariable('Surface_Reflectance', np.float32,
                                               ('size', 'lat', 'lon'))
    group.variables['Surface_Reflectance'][:] = SRF
    Surface_Reflectance.time = date
    Surface_Reflectance.band_names = ['1', '2', '3a']


def SRF_GK2A_nc(SRF, lat, lon, Ins_Dir, name, date):
    Output_path = Ins_Dir + '/' + 'SRF' + name + date + '.nc'
    nc_file = nc.Dataset(Output_path, 'w', format='NETCDF4')
    nc_file.createDimension('lat', len(lat))
    nc_file.createDimension('lon', len(lon))
    nc_file.createDimension('size', SRF.shape[0])

    group = nc_file.createGroup('SRF')

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon

    Surface_Reflectance = group.createVariable('Surface_Reflectance', np.float32,
                                               ('size', 'lat', 'lon'))
    group.variables['Surface_Reflectance'][:] = SRF
    Surface_Reflectance.time = date
    Surface_Reflectance.band_names = ['1', '2', '3', '4', '6']


def algorithm_SRF(TOA_500, sa_z, so_z, A, isCloud, lat, lon, Ins_Dir, name, date):
    atmos_ref_all = list()
    for i in range(7):
        ref_whole = [1, 2, 3, 4, 5, 6, 7]
        table_file = './AT_LUTS/MODIS/MODIS_B' + str(ref_whole[i]) + '.nc'
        atmos_ref = atmos(TOA_500[i], sa_z, so_z, A, table_file, isCloud)
        atmos_ref_all.append(atmos_ref)
    atmos_ref_all = np.array(atmos_ref_all).astype(float)
    return atmos_ref_all


def algorithm_SRF_NOAA19(TOA, sa_z, so_z, A, isCloud, lat, lon, Ins_Dir, name, date):
    atmos_ref_all = list()
    for i in range(3):
        ref_whole = ['1', '2', '3a']
        table_file = './AT_LUTS/AVHRR/AVHRR_B' + ref_whole[i] + '.nc'
        atmos_ref = atmos(TOA[i], sa_z, so_z, A, table_file, isCloud)
        atmos_ref_all.append(atmos_ref)
    atmos_ref_all = np.array(atmos_ref_all).astype(float)
    return atmos_ref_all


def algorithm_SRF_GK2A(TOA, sa_z, so_z, A, isCloud, lat, lon, Ins_Dir, name, date):
    atmos_ref_all = list()
    for i in range(5):
        ref_whole = ['1', '2', '3', '4', '6']
        table_file = './AT_LUTS/AMI/AMI_B' + ref_whole[i] + '.nc'
        atmos_ref = atmos(TOA[i], sa_z, so_z, A, table_file, isCloud)
        atmos_ref_all.append(atmos_ref)
    atmos_ref_all = np.array(atmos_ref_all).astype(float)
    return atmos_ref_all
