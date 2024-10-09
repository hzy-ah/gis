import inspect
import re
import xml.dom.minidom as minidom

import netCDF4 as nc
import numpy as np
from osgeo import gdal


def varname(p):
    for line in inspect.getframeinfo(inspect.currentframe().f_back)[3]:
        m = re.search(r'\bvarname\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)', line)
    if m:
        return m.group(1)


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


def expand_band(data):
    x = data.shape[0]
    y = data.shape[1]
    datax = np.zeros((2 * x, 2 * y))
    for i in range(datax.shape[0]):
        for j in range(datax.shape[1]):
            datax[i][j] = data[int(i / 2)][int(j / 2)]

    return datax


def create_tif(data, out_dir, name):
    data_tif = gdal.GetDriverByName("GTiff")

    data_filename = str(out_dir) + '/' + str(name) + '.tif'

    data_dataset = data_tif.Create(data_filename, data.shape[1], data.shape[0], 1, gdal.GDT_Float32)

    band0 = data_dataset.GetRasterBand(1)

    band0.WriteArray(data)

    return data_filename


def xml_file(data, lat, lon, xml, x, y):
    dom = minidom.getDOMImplementation().createDocument(None, 'VRTDataset', None)
    root = dom.documentElement
    root.setAttribute('rasterXSize', str(y))
    root.setAttribute('rasterYSize', str(x))

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


def fy3c_nc(ref, emissive, cloud, lat_nc, lon_nc, out_dir, year, month, day, out_head1):
    out_path = out_dir + '/' + out_head1 + '_' + str(year) + '_' + str(month) + '_' + str(day) + '.nc'
    f_w = nc.Dataset(out_path, 'w', format='NETCDF4')

    f_w.createDimension('lat', lat_nc.shape[0])
    f_w.createDimension('lon', lon_nc.shape[0])
    f_w.createDimension('size1', ref.shape[0])
    f_w.createDimension('size2', emissive.shape[0])

    group0 = f_w.createGroup('isCloud')
    cloud0 = group0.createVariable('isCloud', np.float32, ('lat', 'lon'))
    group0.variables['isCloud'][:] = cloud
    cloud0.time = str(year) + '_' + str(month) + '_' + str(day)
    group0.createVariable('lat', np.float32, ('lat'))
    group0.variables['lat'][:] = lat_nc
    group0.createVariable('lon', np.float32, ('lon'))
    group0.variables['lon'][:] = lon_nc

    group = f_w.createGroup('TOA')
    reflect = group.createVariable('Top_of_Atmosphere_Reflectance', np.float32, ('size1', 'lat', 'lon'))
    group.variables['Top_of_Atmosphere_Reflectance'][:] = ref
    reflect.time = str(year) + '_' + str(month) + '_' + str(day)
    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat_nc
    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon_nc

    group2 = f_w.createGroup('BT')
    emi = group2.createVariable('Bright_Temperature', np.float32, ('size2', 'lat', 'lon'))
    group2.variables['Bright_Temperature'][:] = emissive
    emi.time = str(year) + '_' + str(month) + '_' + str(day)
    group2.createVariable('lat', np.float32, ('lat'))
    group2.variables['lat'][:] = lat_nc
    group2.createVariable('lon', np.float32, ('lon'))
    group2.variables['lon'][:] = lon_nc


def fy3d_nc(ref, ref2, emissive, cloud, lat_nc, lon_nc,
            lat_nc2, lon_nc2, out_dir, year, month, day, out_head1, sa_z_500, so_z_500, A_500):
    out_path = out_dir + '/' + out_head1 + '_' + str(year) + '_' + str(month) + '_' + str(day) + '.nc'
    f_w = nc.Dataset(out_path, 'w', format='NETCDF4')
    f_w.createDimension('lat_1000', lat_nc.shape[0])
    f_w.createDimension('lon_1000', lon_nc.shape[0])
    f_w.createDimension('lat_500', lat_nc2.shape[0])
    f_w.createDimension('lon_500', lon_nc2.shape[0])
    f_w.createDimension('size1', ref.shape[0])
    f_w.createDimension('size2', ref2.shape[0])
    f_w.createDimension('size3', emissive.shape[0])

    group = f_w.createGroup('TOA_1000')
    reflect = group.createVariable('Top_of_Atmosphere_Reflectance_1000', np.float32, ('size1', 'lat_1000', 'lon_1000'))
    group.variables['Top_of_Atmosphere_Reflectance_1000'][:] = ref
    reflect.time = str(year) + '_' + str(month) + '_' + str(day)
    reflect.band = '5,6,7'
    group.createVariable('lat_1000', np.float32, ('lat_1000'))
    group.variables['lat_1000'][:] = lat_nc

    group.createVariable('lon_1000', np.float32, ('lon_1000'))
    group.variables['lon_1000'][:] = lon_nc

    group2 = f_w.createGroup('TOA_500')
    reflect2 = group2.createVariable('Top_of_Atmosphere_Reflectance_500', np.float32, ('size2', 'lat_500', 'lon_500'))
    group2.variables['Top_of_Atmosphere_Reflectance_500'][:] = ref2
    reflect2.time = str(year) + '_' + str(month) + '_' + str(day)
    reflect2.band = '1,2,3,4'

    group2.createVariable('lat_500', np.float32, ('lat_500'))
    group2.variables['lat_500'][:] = lat_nc2

    group2.createVariable('lon_500', np.float32, ('lon_500'))
    group2.variables['lon_500'][:] = lon_nc2

    group3 = f_w.createGroup('BT_500')
    emi = group3.createVariable('Bright_Temperature_500', np.float32, ('size3', 'lat_500', 'lon_500'))
    group3.variables['Bright_Temperature_500'][:] = emissive
    emi.time = str(year) + '_' + str(month) + '_' + str(day)
    emi.band = '24,25'

    group3.createVariable('lat_500', np.float32, ('lat_500'))
    group3.variables['lat_500'][:] = lat_nc2

    group3.createVariable('lon_500', np.float32, ('lon_500'))
    group3.variables['lon_500'][:] = lon_nc2

    group4 = f_w.createGroup('isCloud')
    iscloud = group4.createVariable('isCloud', np.float32, ('lat_500', 'lon_500'))
    group4.variables['isCloud'][:] = cloud
    iscloud.time = str(year) + '_' + str(month) + '_' + str(day)

    group4.createVariable('lat_500', np.float32, ('lat_500'))
    group4.variables['lat_500'][:] = lat_nc2

    group4.createVariable('lon_500', np.float32, ('lon_500'))
    group4.variables['lon_500'][:] = lon_nc2

    group5 = f_w.createGroup('angle')
    group5.createVariable('SenorZenith', np.float32, ('lat_500', 'lon_500'))
    group5.variables['SenorZenith'][:] = sa_z_500

    group5.createVariable('SolarZenith', np.float32, ('lat_500', 'lon_500'))
    group5.variables['SolarZenith'][:] = so_z_500

    group5.createVariable('RelativeAzimuth', np.float32, ('lat_500', 'lon_500'))
    group5.variables['RelativeAzimuth'][:] = A_500


def jpss_nc(ref, emissive, ref1, emissive1, cloud, lat_nc, lon_nc,
            lat_nc2, lon_nc2, out_dir, year, month, day, out_head1, sa_z, so_z, A):
    out_path = out_dir + '/' + out_head1 + '_' + str(year) + '_' + str(month) + '_' + str(day) + '.nc'
    f_w = nc.Dataset(out_path, 'w', format='NETCDF4')
    f_w.createDimension('lat', lat_nc.shape[0])
    f_w.createDimension('lon', lon_nc.shape[0])
    f_w.createDimension('lat2', lat_nc2.shape[0])
    f_w.createDimension('lon2', lon_nc2.shape[0])
    f_w.createDimension('size1', ref.shape[0])
    f_w.createDimension('size2', emissive.shape[0])
    f_w.createDimension('size3', ref1.shape[0])
    f_w.createDimension('size4', emissive1.shape[0])

    group = f_w.createGroup('TOA_I')
    group.createVariable('Top_of_Atmosphere_Reflectance', np.float32, ('size1', 'lat', 'lon'))
    group.variables['Top_of_Atmosphere_Reflectance'][:] = ref
    group.createVariable('lat_500', np.float32, ('lat'))
    group.variables['lat_500'][:] = lat_nc
    group.createVariable('lon_500', np.float32, ('lon'))
    group.variables['lon_500'][:] = lon_nc

    group2 = f_w.createGroup('BT_I')
    group2.createVariable('Bright_Temperature', np.float32, ('size2', 'lat', 'lon'))
    group2.variables['Bright_Temperature'][:] = emissive
    group2.createVariable('lat_500', np.float32, ('lat'))
    group2.variables['lat_500'][:] = lat_nc
    group2.createVariable('lon_500', np.float32, ('lon'))
    group2.variables['lon_500'][:] = lon_nc

    group3 = f_w.createGroup('TOA_M')
    group3.createVariable('Top_of_Atmosphere_Reflectance', np.float32, ('size3', 'lat2', 'lon2'))
    group3.variables['Top_of_Atmosphere_Reflectance'][:] = ref1
    group3.createVariable('lat_1000', np.float32, ('lat2'))
    group3.variables['lat_1000'][:] = lat_nc2
    group3.createVariable('lon_1000', np.float32, ('lon2'))
    group3.variables['lon_1000'][:] = lon_nc2

    group4 = f_w.createGroup('BT_M')
    group4.createVariable('Bright_Temperature', np.float32, ('size4', 'lat2', 'lon2'))
    group4.variables['Bright_Temperature'][:] = emissive1
    group4.createVariable('lat_1000', np.float32, ('lat2'))
    group4.variables['lat_1000'][:] = lat_nc2
    group4.createVariable('lon_1000', np.float32, ('lon2'))
    group4.variables['lon_1000'][:] = lon_nc2

    group5 = f_w.createGroup('isCloud')
    group5.createVariable('isCloud', np.float32, ('lat', 'lon'))
    group5.variables['isCloud'][:] = cloud
    group5.createVariable('lat_500', np.float32, ('lat'))
    group5.variables['lat_500'][:] = lat_nc
    group5.createVariable('lon_500', np.float32, ('lon'))
    group5.variables['lon_500'][:] = lon_nc

    group6 = f_w.createGroup('angle')
    group6.createVariable('SenorZenith', np.float32, ('lat', 'lon'))
    group6.variables['SenorZenith'][:] = sa_z

    group6.createVariable('SolarZenith', np.float32, ('lat', 'lon'))
    group6.variables['SolarZenith'][:] = so_z

    group6.createVariable('RelativeAzimuth', np.float32, ('lat', 'lon'))
    group6.variables['RelativeAzimuth'][:] = A


def atmos_fy3c_nc(srf, lat_nc, lon_nc, out_dir, year, month, day, out_head, text):
    out_path = out_dir + '/' + out_head + '_' + str(year) + '_' + str(month) + '_' + str(day) + '.nc'
    f_w = nc.Dataset(out_path, 'w', format='NETCDF4')
    f_w.createDimension('lat', lat_nc.shape[0])
    f_w.createDimension('lon', lon_nc.shape[0])
    f_w.createDimension('size', srf.shape[0])

    group = f_w.createGroup('SRF')
    data = group.createVariable('Surface_Reflectance', np.float32, ('size', 'lat', 'lon'))
    group.variables['Surface_Reflectance'][:] = srf
    data.band = text
    data.time = str(year) + '_' + str(month) + '_' + str(day)

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat_nc

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon_nc

    f_w.close()


def atmos_fy3d_nc(ref_1000, ref_500, lat_nc_1000, lon_nc_1000,
                  lat_nc_500, lon_nc_500, out_dir, year, month, day, out_head1, sa_z_500, so_z_500, A_500):
    out_path = out_dir + '/' + out_head1 + '_' + str(year) + '_' + str(month) + '_' + str(day) + '.nc'
    f_w = nc.Dataset(out_path, 'w', format='NETCDF4')
    f_w.createDimension('lat_1000', lat_nc_1000.shape[0])
    f_w.createDimension('lon_1000', lon_nc_1000.shape[0])
    f_w.createDimension('lat_500', lat_nc_500.shape[0])
    f_w.createDimension('lon_500', lon_nc_500.shape[0])
    f_w.createDimension('size_1000', ref_1000.shape[0])
    f_w.createDimension('size_500', ref_500.shape[0])

    group = f_w.createGroup('SRF_1000')
    reflect_1000 = group.createVariable('Surface_Reflectance_1000', np.float32, ('size_1000', 'lat_1000', 'lon_1000'))
    group.variables['Surface_Reflectance_1000'][:] = ref_1000
    reflect_1000.time = str(year) + '_' + str(month) + '_' + str(day)
    reflect_1000.band = '6,7'

    group.createVariable('lat_1000', np.float32, ('lat_1000'))
    group.variables['lat_1000'][:] = lat_nc_1000

    group.createVariable('lon_1000', np.float32, ('lon_1000'))
    group.variables['lon_1000'][:] = lon_nc_1000

    group2 = f_w.createGroup('SRF_500')
    reflect_500 = group2.createVariable('Surface_Reflectance_500', np.float32, ('size_500', 'lat_500', 'lon_500'))
    group2.variables['Surface_Reflectance_500'][:] = ref_500
    reflect_500.time = str(year) + '_' + str(month) + '_' + str(day)
    reflect_500.band = '1,2,3,4'

    group2.createVariable('lat_500', np.float32, ('lat_500'))
    group2.variables['lat_500'][:] = lat_nc_500

    group2.createVariable('lon_500', np.float32, ('lon_500'))
    group2.variables['lon_500'][:] = lon_nc_500

    group5 = f_w.createGroup('angle')
    group5.createVariable('SenorZenith', np.float32, ('lat_500', 'lon_500'))
    group5.variables['SenorZenith'][:] = sa_z_500

    group5.createVariable('SolarZenith', np.float32, ('lat_500', 'lon_500'))
    group5.variables['SolarZenith'][:] = so_z_500

    group5.createVariable('RelativeAzimuth', np.float32, ('lat_500', 'lon_500'))
    group5.variables['RelativeAzimuth'][:] = A_500


def atmos_jpss_nc(srf, srf2, lat_nc, lon_nc, lat_nc2, lon_nc2,
                  out_dir, year, month, day, out_head, text, sa_z, so_z, A):
    out_path = out_dir + '/' + out_head + '_' + str(year) + '_' + str(month) + '_' + str(day) + '.nc'
    f_w = nc.Dataset(out_path, 'w', format='NETCDF4')
    f_w.createDimension('lat', lat_nc.shape[0])
    f_w.createDimension('lon', lon_nc.shape[0])
    f_w.createDimension('lat2', lat_nc2.shape[0])
    f_w.createDimension('lon2', lon_nc2.shape[0])
    f_w.createDimension('size', srf.shape[0])
    f_w.createDimension('size2', srf2.shape[0])

    group = f_w.createGroup('SRF_I')
    data = group.createVariable('Surface_Reflectance', np.float32, ('size', 'lat', 'lon'))
    group.variables['Surface_Reflectance'][:] = srf

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat_nc

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon_nc

    data.band = text[0]
    data.time = str(year) + '_' + str(month) + '_' + str(day)

    group2 = f_w.createGroup('SRF_M')
    data2 = group2.createVariable('Surface_Reflectance', np.float32, ('size2', 'lat2', 'lon2'))
    group2.variables['Surface_Reflectance'][:] = srf2

    group2.createVariable('lat', np.float32, ('lat'))
    group2.variables['lat'][:] = lat_nc

    group2.createVariable('lon', np.float32, ('lon'))
    group2.variables['lon'][:] = lon_nc

    data2.band = text[1]
    data2.time = str(year) + '_' + str(month) + '_' + str(day)

    group6 = f_w.createGroup('angle')
    group6.createVariable('SenorZenith', np.float32, ('lat', 'lon'))
    group6.variables['SenorZenith'][:] = sa_z

    group6.createVariable('SolarZenith', np.float32, ('lat', 'lon'))
    group6.variables['SolarZenith'][:] = so_z

    group6.createVariable('RelativeAzimuth', np.float32, ('lat', 'lon'))
    group6.variables['RelativeAzimuth'][:] = A

    f_w.close()


def zhishu_nc(zhushu, lat_nc, lon_nc, out_dir, year, month, day, out_head, name):
    out_path = out_dir + '/' + name + out_head + '_' + str(year) + '_' + str(month) + '_' + str(day) + '.nc'
    f_w = nc.Dataset(out_path, 'w', format='NETCDF4')
    f_w.createDimension('lat', lat_nc.shape[0])
    f_w.createDimension('lon', lon_nc.shape[0])

    group = f_w.createGroup(name)

    data = group.createVariable(name, np.float32, ('lat', 'lon'))
    group.variables[name][:] = zhushu
    data.time = str(year) + '_' + str(month) + '_' + str(day)

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat_nc

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon_nc

    f_w.close()


def zs_month(data, out_dir, header, lat, lon):
    data_aver = -999 * np.ones((data.shape[1], data.shape[2]))
    for i in range(data.shape[1]):
        for j in range(data.shape[2]):
            dd = []
            for num in range(data.shape[0]):
                if data[num][i][j] != -999:
                    dd.append(data[num][i][j])

            if dd:
                data_aver[i][j] = np.median(np.array(dd))

    out_path = out_dir + '/' + header[0] + '_' + header[2] + '.nc'
    f_w = nc.Dataset(out_path, 'w', format='NETCDF4')
    f_w.createDimension('1', data_aver.shape[0])
    f_w.createDimension('2', data_aver.shape[1])
    f_w.createDimension('lat', len(lat))
    f_w.createDimension('lon', len(lon))

    group = f_w.createGroup(header[1])
    data = group.createVariable(header[1], np.float32, ('1', '2'))
    group.variables[header[1]][:] = data_aver

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon

    return data_aver


def zs_year(data, out_dir, header, lat, lon):
    data_aver = np.zeros((data.shape[1], data.shape[2]))
    for i in range(data.shape[1]):
        for j in range(data.shape[2]):
            dd = []
            for num in range(data.shape[0]):

                if data[num][i][j] != -999:
                    dd.append(data[num][i][j])

            if dd:
                data_aver[i][j] = np.median(np.array(dd))

    out_path = out_dir + '/' + header[0] + '_' + header[2] + '.nc'

    f_w = nc.Dataset(out_path, 'w', format='NETCDF4')
    f_w.createDimension('1', data_aver.shape[0])
    f_w.createDimension('2', data_aver.shape[1])
    f_w.createDimension('lat', len(lat))
    f_w.createDimension('lon', len(lon))

    group = f_w.createGroup(header[1])
    data = group.createVariable(header[1], np.float32, ('1', '2'))
    group.variables[header[1]][:] = data_aver

    group.createVariable('lat', np.float32, ('lat'))
    group.variables['lat'][:] = lat

    group.createVariable('lon', np.float32, ('lon'))
    group.variables['lon'][:] = lon

    return data_aver


def replace_fill_value(image):
    neighbors = [(0, -1), (-1, 0), (0, 1), (1, 0), (1, 1), (1, -1), (-1, -1), (-1, 1)]
    image2 = image.copy()
    for i in range(5, image.shape[0] - 5):
        for j in range(5, image.shape[1] - 5):
            if image[i][j] == -999:
                valid_pixels = []
                for dx, dy in neighbors:
                    x = i + dx
                    y = j + dy

                    if x >= 0 and x < image.shape[0] and y >= 0 and y < image.shape[1]:
                        pixel = image[x, y]

                        if pixel != -999:
                            valid_pixels.append(pixel)

                if valid_pixels:
                    average = np.mean(valid_pixels)
                    image2[i][j] = average

    return image2


def remove_rows_with_value(matrix, value):
    # 创建一个布尔索引，表示每一行是否全是给定的value
    row_mask = np.all(matrix < value, axis=1)

    # 使用布尔索引删除满足条件的行
    new_matrix = matrix[~row_mask]

    # 返回删除的行的索引位置
    removed_rows = np.where(row_mask)[0]

    return new_matrix, removed_rows
