import numpy as np
import netCDF4 as nc


def GetDataList_AVHRR_AMI(file_list_all):
    file_list_RVI = []
    file_list_NDVI = []
    file_list_SAVI = []
    file_list_EVI = []
    file_list_LAI = []

    for file in file_list_all:
        if file.find('RVI') != -1:
            file_list_RVI.append(file)
        elif file.find('NDVI') != -1:
            file_list_NDVI.append(file)
        elif file.find('SAVI') != -1:
            file_list_SAVI.append(file)
        elif file.find('EVI') != -1:
            file_list_EVI.append(file)
        elif file.find('LAI') != -1:
            file_list_LAI.append(file)

    return file_list_RVI, file_list_NDVI, file_list_SAVI, file_list_EVI, file_list_LAI


def GetDataList(file_list_all):
    file_list_RVI = []
    file_list_NDVI = []
    file_list_SAVI = []
    file_list_GVI = []
    file_list_EVI = []
    file_list_LAI = []

    for file in file_list_all:
        if file.find('RVI') != -1:
            file_list_RVI.append(file)
        elif file.find('NDVI') != -1:
            file_list_NDVI.append(file)
        elif file.find('SAVI') != -1:
            file_list_SAVI.append(file)
        elif file.find('GVI') != -1:
            file_list_GVI.append(file)
        elif file.find('EVI') != -1:
            file_list_EVI.append(file)
        elif file.find('LAI') != -1:
            file_list_LAI.append(file)

    return file_list_RVI, file_list_NDVI, file_list_SAVI, file_list_GVI, file_list_EVI, file_list_LAI


def CalData(Ins_Dir, file_list, Type, size):
    data_mon = np.ndarray(shape=[len(file_list), size, size])
    lat = nc.Dataset(Ins_Dir + file_list[0]).groups[Type].variables['lat'][:]
    lon = nc.Dataset(Ins_Dir + file_list[0]).groups[Type].variables['lon'][:]

    for i in range(len(file_list)):
        NC = nc.Dataset(Ins_Dir + file_list[i])
        group = NC.groups[Type]
        data_mon[i] = group.variables[Type][:]

    data = data_mon[0]
    for i in range(data_mon.shape[1]):
        for j in range(data_mon.shape[2]):
            data_list = []
            for k in range(data_mon.shape[0]):
                if data[k][i][j] != -999:
                    data_list.append(data[k][i][j])
            data[i][j] = np.mean(data_list)

    return data, lat, lon


def Write_nc(data, lat, lon, Mon_Dir, Type, name, date):
    Output_path = Mon_Dir + '/' + 'Mon-' + Type + name + date + '.nc'
    nc_file = nc.Dataset(Output_path, 'w', format='NETCDF4')
    nc_file.createDimension('lat', len(lat))
    nc_file.createDimension('lon', len(lon))
    nc_file.createDimension('size', 1)

    group = nc_file.createGroup(Type)

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon

    dataset = group.createVariable(Type, np.float32,
                                   ('size', 'lat', 'lon'))
    group.variables[Type][:] = data
    dataset.time = date
    dataset.band_names = [Type]


def Write_nc_Day(data, lat, lon, Mon_Dir, Type, name, date):
    Output_path = Mon_Dir + '/' + 'Day-' + Type + name + date + '.nc'
    nc_file = nc.Dataset(Output_path, 'w', format='NETCDF4')
    nc_file.createDimension('lat', len(lat))
    nc_file.createDimension('lon', len(lon))
    nc_file.createDimension('size', 1)

    group = nc_file.createGroup(Type)

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon

    dataset = group.createVariable(Type, np.float32,
                                   ('size', 'lat', 'lon'))
    group.variables[Type][:] = data
    dataset.time = date
    dataset.band_names = [Type]
