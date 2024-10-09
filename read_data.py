import numpy as np
import math
import os
import struct
import create_data as cd
import sys
from tqdm import tqdm


def read_geo(geo, data_type):
    # 具体的波段读取
    if data_type == 'VIIRS':
        allkey = []
        for key in geo.keys():
            allkey.append(key)

        allkey2 = []
        for key in geo[allkey[0]].keys():
            allkey2.append(key)

        h5_geo_h = geo[allkey[0]][allkey2[0]]['Height']
        h5_geo_Lat = geo[allkey[0]][allkey2[0]]['Latitude']
        h5_geo_Lon = geo[allkey[0]][allkey2[0]]['Longitude']
        SatelliteAzimuthAngle = geo[allkey[0]][allkey2[0]]['SatelliteAzimuthAngle']
        SatelliteZenithAngle = geo[allkey[0]][allkey2[0]]['SatelliteZenithAngle']
        SolarAzimuthAngle = geo[allkey[0]][allkey2[0]]['SolarAzimuthAngle']
        SolarZenithAngle = geo[allkey[0]][allkey2[0]]['SolarZenithAngle']

        h = np.array(h5_geo_h).astype(float)
        lat = np.array(h5_geo_Lat).astype(float)
        lon = np.array(h5_geo_Lon).astype(float)
        sa_a = np.array(SatelliteAzimuthAngle).astype(float)
        sa_z = np.array(SatelliteZenithAngle).astype(float)
        so_a = np.array(SolarAzimuthAngle).astype(float)
        so_z = np.array(SolarZenithAngle).astype(float)

        return h, lat, lon, sa_a, sa_z, so_a, so_z
    if data_type == 'FY3C':
        allkey = []
        for key in geo.keys():
            allkey.append(key)
        h5_geo_h = geo[allkey[0]]['DEM']
        h5_geo_Lat = geo[allkey[0]]['Latitude']
        h5_geo_Lon = geo[allkey[0]]['Longitude']
        SatelliteAzimuthAngle = geo[allkey[0]]['SensorAzimuth']
        SatelliteZenithAngle = geo[allkey[0]]['SensorZenith']
        SolarAzimuthAngle = geo[allkey[0]]['SolarAzimuth']
        SolarZenithAngle = geo[allkey[0]]['SolarZenith']
        SatelliteAzimuthAngle = np.array(SatelliteAzimuthAngle).astype(float) * 0.01
        SatelliteZenithAngle = np.array(SatelliteZenithAngle).astype(float) * 0.01
        SolarAzimuthAngle = np.array(SolarAzimuthAngle).astype(float) * 0.01
        SolarZenithAngle = np.array(SolarZenithAngle).astype(float) * 0.01

        h = np.array(h5_geo_h).astype(float)
        lat = np.array(h5_geo_Lat).astype(float)
        lon = np.array(h5_geo_Lon).astype(float)
        sa_a = np.array(SatelliteAzimuthAngle).astype(float)
        sa_z = np.array(SatelliteZenithAngle).astype(float)
        so_a = np.array(SolarAzimuthAngle).astype(float)
        so_z = np.array(SolarZenithAngle).astype(float)

        return h, lat, lon, sa_a, sa_z, so_a, so_z
    if data_type == 'FY3D_1000':
        h5_geo_h = geo['Geolocation']['DEM']
        h5_geo_Lat = geo['Geolocation']['Latitude']
        h5_geo_Lon = geo['Geolocation']['Longitude']
        SatelliteAzimuthAngle = geo['Geolocation']['SensorAzimuth']
        SatelliteZenithAngle = geo['Geolocation']['SensorZenith']
        SolarAzimuthAngle = geo['Geolocation']['SolarAzimuth']
        SolarZenithAngle = geo['Geolocation']['SolarZenith']
        SatelliteAzimuthAngle = np.array(SatelliteAzimuthAngle).astype(float) * 0.01
        SatelliteZenithAngle = np.array(SatelliteZenithAngle).astype(float) * 0.01
        SolarAzimuthAngle = np.array(SolarAzimuthAngle).astype(float) * 0.01
        SolarZenithAngle = np.array(SolarZenithAngle).astype(float) * 0.01

        h = np.array(h5_geo_h).astype(float)
        lat = np.array(h5_geo_Lat).astype(float)
        lon = np.array(h5_geo_Lon).astype(float)
        sa_a = np.array(SatelliteAzimuthAngle).astype(float)
        sa_z = np.array(SatelliteZenithAngle).astype(float)
        so_a = np.array(SolarAzimuthAngle).astype(float)
        so_z = np.array(SolarZenithAngle).astype(float)

        return h, lat, lon, sa_a, sa_z, so_a, so_z
    if data_type == 'FY3D_250':

        h5_geo_Lat = geo['Latitude']
        h5_geo_Lon = geo['Longitude']

        lat = np.array(h5_geo_Lat).astype(float)
        lon = np.array(h5_geo_Lon).astype(float)

        lat = lat[::2, ::2]
        lon = lon[::2, ::2]

        x_max = 0
        x_min = 10000
        y_max = 0
        y_min = 10000

        for i in range(lat.shape[0]):
            for j in range(lat.shape[1]):
                if lat[i][j] >= 20 and lat[i][j] <= 30 and lon[i][j] >= 114 and lon[i][j] <= 124:
                    if i > x_max:
                        x_max = i

                    if i < x_min:
                        x_min = i

                    if j > y_max:
                        y_max = j

                    if j < y_min:
                        y_min = j

        data_range = (x_min, x_max, y_min, y_max)

        lat = lat[x_min:x_max + 1, y_min:y_max + 1]
        lon = lon[x_min:x_max + 1, y_min:y_max + 1]

        return lat, lon, data_range


def read_band(data, data_type, band_type=0, data_range=None):
    if data_type == 'VIIRS':
        if band_type == 'ref_i' or band_type == 'ref_m':
            allkey = []
            for key in data.keys():
                allkey.append(key)
            allkey2 = []
            for key in data[allkey[0]].keys():
                allkey2.append(key)
            Ref = data[allkey[0]][allkey2[0]]['Reflectance']
            factor = data[allkey[0]][allkey2[0]]['ReflectanceFactors']
            Ref = np.array(Ref).astype(float)
            factor = np.array(factor).astype(float)
            if band_type == 'ref_i':
                Ref = viirs_i(Ref)

            if band_type == 'ref_m':
                Ref = viirs_m(Ref)

            if data_range is not None:
                x_min = data_range[0]
                x_max = data_range[1]
                y_min = data_range[2]
                y_max = data_range[3]
                Ref = Ref[x_min:x_max + 1, y_min:y_max + 1]

            for i in range(Ref.shape[0]):
                for j in range(Ref.shape[1]):
                    Ref[i][j] = Ref[i][j] * factor[0] + factor[1]

            return Ref

        if band_type == 'bt_i' or band_type == 'bt_m':
            allkey = []
            for key in data.keys():
                allkey.append(key)
            allkey2 = []
            for key in data[allkey[0]].keys():
                allkey2.append(key)
            Bt = data[allkey[0]][allkey2[0]]['BrightnessTemperature']
            factor = data[allkey[0]][allkey2[0]]['BrightnessTemperatureFactors']
            Bt = np.array(Bt).astype(float)
            factor = np.array(factor).astype(float)

            if band_type == 'bt_i':
                Bt = viirs_i(Bt)

            if band_type == 'bt_m':
                Bt = viirs_m(Bt)
            # # 去除条带
            # tag = np.zeros(Bt.shape[0])
            # tag2 = np.zeros(Bt.shape[0])
            # x = 0
            # y = 0
            # for i in range(10, Bt.shape[0] - 10):
            #     if Bt[i][1] == 65533 and Bt[i][20] == 65533:
            #         tag[i] = 1
            #
            # for j in range(10, tag.shape[0] - 10):
            #     if tag[j] == 1 and tag[j - 1] == 0:
            #         tag2[j] = 1
            #         if x == 0:
            #             x = j
            #     if tag[j] == 1 and tag[j + 1] == 0:
            #         tag2[j] = 2
            #         if y == 0:
            #             y = j
            #
            # size = int(math.fabs(y - x + 1))
            # qsize = int(size / 4)
            #
            # for i in range(5, Bt.shape[0]):
            #     for j in range(Bt.shape[1]):
            #         if tag2[i] == 1:
            #             if Bt[i][j] == 65533:
            #                 for k in range(qsize):
            #                     Bt[i + k][j] = Bt[i - qsize + k][j]
            #                     Bt[i + k + 3 * qsize][j] = Bt[i + size + k][j]
            #             if Bt[i + qsize][j] == 65533:
            #                 for k in range(qsize):
            #                     Bt[i + qsize + k][j] = Bt[i + k][j]
            #                     Bt[i + k + 2 * qsize][j] = Bt[i + 3 * qsize + k][j]

            if data_range is not None:
                x_min = data_range[0]
                x_max = data_range[1]
                y_min = data_range[2]
                y_max = data_range[3]

                Bt = Bt[x_min:x_max + 1, y_min:y_max + 1]

            # 辐射定标
            for i in range(Bt.shape[0]):
                for j in range(Bt.shape[1]):
                    Bt[i][j] = Bt[i][j] * factor[0] + factor[1]

            return Bt

    if data_type == 'FY3C':
        drn = data.attrs['Day Or Night Flag']
        if drn[0] != 68:
            print('Can not find the idea FY3C-VIRR L1B file at your input date and time!')

            sys.exit()

        ref_whole = [1, 2, 6, 7, 8, 9, 10]
        rad_whole = [3, 4, 5]

        Ref = data['Data']['EV_RefSB']
        Coeffcients = data.attrs['RefSB_Cal_Coefficients']
        E0 = data.attrs['RefSB_Solar_Irradiance']

        Ref = np.array(Ref).astype(float)
        Coeffcients = np.array(Coeffcients).astype(float)
        E0 = np.array(E0).astype(float)

        for num in range(Ref.shape[0]):

            k = Coeffcients[2 * num]
            b = Coeffcients[2 * num + 1]
            for i in range(Ref.shape[1]):
                for j in range(Ref.shape[2]):
                    Ref[num][i][j] = (Ref[num][i][j] * k + b) / 100

        Emissive = data['Data']['EV_Emissive']
        Scale = data['Data']['Emissive_Radiance_Scales']
        Offset = data['Data']['Emissive_Radiance_Offsets']
        bn = data.attrs['Prelaunch_Nonlinear_Coefficients']
        wave_n = data.attrs['Emissive_Centroid_Wave_Number']
        bt = data.attrs['Emissive_BT_Coefficients']

        Emissive = np.array(Emissive).astype(float)
        Scale = np.array(Scale).astype(float)
        Offset = np.array(Offset).astype(float)

        bn = np.array(bn).astype(float)
        wave_n = np.array(wave_n).astype(float)
        bt = np.array(bt).astype(float)
        for num in range(Emissive.shape[0]):

            b0 = bn[3 * num]
            b1 = bn[3 * num + 1]
            b2 = bn[3 * num + 2]

            c1 = 1.1910427 * math.pow(10, -5)
            c2 = 1.4387752
            vc = wave_n[num]

            A = bt[2 * num]
            B = bt[2 * num + 1]

            for i in range(Emissive.shape[1]):
                for j in range(Emissive.shape[2]):
                    Emissive[num][i][j] = Scale[i][num] * Emissive[num][i][j] + Offset[i][num]
                    Emissive[num][i][j] = b0 + (1 + b1) * Emissive[num][i][j] + b2 * Emissive[num][i][j] * \
                                          Emissive[num][i][j]
                    if Emissive[num][i][j] < 0:
                        Emissive[num][i][j] = 280

                    else:
                        Emissive[num][i][j] = float(c2 * vc) / (
                            math.log(1 + (float(c1 * math.pow(vc, 3)) / float(Emissive[num][i][j]))))
                        Emissive[num][i][j] = (Emissive[num][i][j] - A) / B

        return Ref, Emissive

    if data_type == 'FY3D_1000':

        drn = data.attrs['Day Or Night Flag']

        if drn[0] != 68:
            print('Can not find the idea observation at your input date and time!')
            sys.exit()

        ref_whole = [5, 6, 7]

        Ref = data['Data']['EV_1KM_RefSB']
        Coeffcients = data['Calibration']['VIS_Cal_Coeff']

        Ref = np.array(Ref).astype(float)[:3]

        for num in range(len(ref_whole)):
            # print('计算第' + str(ref_whole[num]) + '波段')
            cal0 = Coeffcients[num + 4][0]
            cal1 = Coeffcients[num + 4][1]
            cal2 = Coeffcients[num + 4][2]
            for i in range(Ref.shape[1]):
                for j in range(Ref.shape[2]):
                    Ref[num][i][j] = (cal2 * Ref[num][i][j] * Ref[num][i][j] + cal1 * Ref[num][i][j] + cal0) / 100

        return Ref

    if data_type == 'FY3D_250':
        ref1 = data['Data']['EV_250_RefSB_b1']
        ref2 = data['Data']['EV_250_RefSB_b2']
        ref3 = data['Data']['EV_250_RefSB_b3']
        ref4 = data['Data']['EV_250_RefSB_b4']
        Coeffcients = data['Calibration']['VIS_Cal_Coeff']

        ref_whole = [1, 2, 3, 4]
        Ref = list()
        Ref.append(ref1)
        Ref.append(ref2)
        Ref.append(ref3)
        Ref.append(ref4)

        Ref = np.array(Ref).astype(float)
        Ref = Ref[:, ::2, ::2]

        x_min = data_range[0]
        x_max = data_range[1]
        y_min = data_range[2]
        y_max = data_range[3]

        data_range = [x_min, x_max, y_min, y_max]
        Ref = Ref[:, x_min:x_max + 1, y_min:y_max + 1]

        for num in range(len(ref_whole)):

            cal0 = Coeffcients[num][0]
            cal1 = Coeffcients[num][1]
            cal2 = Coeffcients[num][2]
            for i in range(Ref.shape[1]):
                for j in range(Ref.shape[2]):
                    Ref[num][i][j] = (cal2 * Ref[num][i][j] * Ref[num][i][j] + cal1 * Ref[num][i][j] + cal0) / 100

        emissive24 = data['Data']['EV_250_Emissive_b24']
        emissive25 = data['Data']['EV_250_Emissive_b25']
        A = data.attrs['TBB_Trans_Coefficient_A']
        B = data.attrs['TBB_Trans_Coefficient_B']
        slope = emissive24.attrs['Slope']
        intercept = emissive24.attrs['Intercept']
        slope2 = emissive25.attrs['Slope']
        intercept2 = emissive25.attrs['Intercept']

        wave_n = [2634.359, 2471.654, 1382.621, 1168.182, 933.364, 836.941]
        rad_whole = [24, 25]
        rad_index = [4, 5]

        Emissive = list()
        Emissive.append(emissive24)
        Emissive.append(emissive25)
        Emissive = np.array(Emissive).astype(float)
        Slope = list()
        Intercept = list()
        Slope.append(slope)
        Slope.append(slope2)
        Intercept.append(intercept)
        Intercept.append(intercept2)

        Emissive = Emissive[:, ::2, ::2]
        Emissive = Emissive[:, x_min:x_max + 1, y_min:y_max + 1]

        for num in range(Emissive.shape[0]):

            a = A[rad_index[num]]
            b = B[rad_index[num]]

            c1 = 1.1910427 * math.pow(10, -5)
            c2 = 1.4387752
            vc = wave_n[rad_index[num]]
            for i in range(Emissive.shape[1]):
                for j in range(Emissive.shape[2]):
                    if Emissive[num][i][j] <= 0:
                        Emissive[num][i][j] = 280
                    else:
                        Emissive[num][i][j] = Emissive[num][i][j] * Slope[num] + Intercept[num]
                        Emissive[num][i][j] = c2 * vc / (math.log(1 + (c1 * math.pow(vc, 3) / Emissive[num][i][j])))
                        Emissive[num][i][j] = Emissive[num][i][j] * a + b

        return Ref, Emissive, data_range


def data_half(data):
    dataf = np.zeros((data.shape[0], int(data.shape[1] / 2), int(data.shape[2] / 2)))
    for num in range(dataf.shape[0]):
        for i in range(dataf.shape[1]):
            for j in range(dataf.shape[2]):
                dataf[num][i][j] = data[num][i * 2][j * 2]

    return dataf


def data_half2(data):
    dataf = np.zeros((int(data.shape[0] / 2), int(data.shape[1] / 2)))
    for i in range(dataf.shape[0]):
        for j in range(dataf.shape[1]):
            dataf[i][j] = data[i * 2][j * 2]
    return dataf


def relative_angle(sa_a, sa_z, so_a, so_z):
    dtr = math.pi / 180
    rtd = 180 / math.pi

    A = np.zeros((sa_z.shape[0], sa_z.shape[1]))
    for i in range(sa_z.shape[0]):
        for j in range(sa_z.shape[1]):
            raa = math.fabs(sa_a[i][j] - so_a[i][j])
            if raa > 180:
                raa = 360 - raa
            A[i][j] = math.acos(math.cos(sa_z[i][j] * dtr) * math.cos(so_z[i][j] * dtr) +
                                math.sin(sa_z[i][j] * dtr) * math.sin(so_z[i][j] * dtr) *
                                math.cos(raa * dtr)) * rtd
    return A


def read_lat_lon(ds):
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    geotransform = ds.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    lat = np.zeros(rows)
    lon = np.zeros(cols)

    for i in range(rows):
        lat[i] = originY + pixelHeight * i

    for i in range(cols):
        lon[i] = originX + pixelWidth * i

    return lat, lon


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

    ir1[0] = ir1[0] / 10 ** 8
    ir1[1] = ir1[1] / 10 ** 6
    ir1[2] = ir1[2] / 10 ** 6

    ir2[0] = ir2[0] / 10 ** 8
    ir2[1] = ir2[1] / 10 ** 6
    ir2[2] = ir2[2] / 10 ** 6

    ir3[0] = ir3[0] / 10 ** 8
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

    for num in tqdm(range(int(size / 22016))):

        if num == 0:
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

            ds.seek(num * skip + 228, 0)
            data = ds.read(72)
            ir_all = list()
            for i in range(18):
                ir_s = (struct.unpack('i', data[4 * i:4 * i + 4])[0])
                ir_all.append(ir_s)

            ir = noaa_ir(ir_all)

            ds.seek(num * skip + 328, 0)
            data = ds.read(306)
            angle_all = list()
            for i in range(153):
                angle_s = float(struct.unpack('h', data[2 * i:2 * i + 2])[0])
                angle_all.append(angle_s)

            so, sa, a = noaa_angle(angle_all)

            so_expand = cd.expand_data(so)
            sa_expand = cd.expand_data(sa)
            a_expand = cd.expand_data(a)

            so_all.append(so_expand)
            sa_all.append(sa_expand)
            a_all.append(a_expand)

            ds.seek(num * skip + 640, 0)
            data = ds.read(408)
            geo_all = list()
            for i in range(102):
                geo_s = float(struct.unpack('i', data[4 * i:4 * i + 4])[0])
                geo_all.append(geo_s)

            lat, lon = geo_avhrr(geo_all)

            lat_expand = cd.expand_data(lat)
            lon_expand = cd.expand_data(lon)

            lat_all.append(lat_expand)
            lon_all.append(lon_expand)

            ds.seek(num * skip + 1264, 0)
            data = ds.read(20480)
            band_all = list()
            for i in range(10240):
                band_s = float(struct.unpack('h', data[2 * i:2 * i + 2])[0])
                band_all.append(band_s)

            band1, band2, band3, band4, band5 = band_avhrr(band_all)

            c1 = 1.1910427 / 10 ** 5
            c2 = 1.4387752

            if num == 2:
                print(band1)

            for i in range(band1.shape[0]):
                if band1[i] > col[0][4]:
                    band1[i] = band1[i] * col[0][2] + col[0][3]


                else:
                    band1[i] = band1[i] / 10
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
            if num == 2:
                print(band1)

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
    ref.append(band1_all)
    ref.append(band2_all)
    ref.append(band3_all)

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

    return ref, irt, sa_all, so_all, a_all, lat_all, lon_all


def viirs_i(Ref):
    tow_tie_sign = 0
    short_tie_length = 0
    long_tie_length = 0
    short_tie_length2 = 0
    long_tie_length2 = 0
    for i in range(Ref.shape[0]):
        if i > 4:
            if Ref[i][2] == 65533:
                for j in range(Ref.shape[1] - 100):
                    if Ref[i][j + 100] != 65533:
                        short_tie_length = j + 100
                        break
                for j in range(Ref.shape[1] - 100):
                    if Ref[i][-j - 100] != 65533:
                        short_tie_length2 = j + 100
                        break
                for j in range(Ref.shape[1] - 100):
                    if Ref[i + 2][j + 100] != 65533:
                        long_tie_length = j + 100
                        break

                for j in range(Ref.shape[1] - 100):
                    if Ref[i + 2][-j - 100] != 65533:
                        long_tie_length2 = j + 100
                        break

                tow_tie_sign = 1

        if tow_tie_sign == 1:
            break

    for i in range(Ref.shape[0]):
        if i > 4 and i < (Ref.shape[0] - 4):
            if Ref[i][0] == 65533 and Ref[i - 1][0] != 65533:
                for j in range(short_tie_length):
                    Ref[i][j] = Ref[i - 1][j]
                    Ref[i + 1][j] = Ref[i - 1][j]
                    Ref[i + 6][j] = Ref[i + 8][j]
                    Ref[i + 7][j] = Ref[i + 8][j]

                for j in range(short_tie_length2):
                    Ref[i][-j] = Ref[i - 1][-j]
                    Ref[i + 1][-j] = Ref[i - 1][-j]
                    Ref[i + 6][-j] = Ref[i + 8][-j]
                    Ref[i + 7][-j] = Ref[i + 8][-j]

                for j in range(long_tie_length):
                    Ref[i + 2][j] = Ref[i][j]
                    Ref[i + 3][j] = Ref[i][j]
                    Ref[i + 4][j] = Ref[i + 8][j]
                    Ref[i + 5][j] = Ref[i + 8][j]

                for j in range(long_tie_length2):
                    Ref[i + 2][-j] = Ref[i][-j]
                    Ref[i + 3][-j] = Ref[i][-j]
                    Ref[i + 4][-j] = Ref[i + 8][-j]
                    Ref[i + 5][-j] = Ref[i + 8][-j]

    return Ref


def viirs_m(Ref):
    tow_tie_sign = 0
    short_tie_length = 0
    long_tie_length = 0
    short_tie_length2 = 0
    long_tie_length2 = 0
    for i in range(Ref.shape[0]):
        if i > 4:
            if Ref[i][2] == 65533:

                for j in range(Ref.shape[1] - 100):
                    if Ref[i][j + 100] != 65533:
                        short_tie_length = j + 100
                        break
                for j in range(Ref.shape[1] - 100):
                    if Ref[i][-j - 100] != 65533:
                        short_tie_length2 = j + 100
                        break
                for j in range(Ref.shape[1] - 100):
                    if Ref[i + 1][j + 100] != 65533:
                        long_tie_length = j + 100
                        break

                for j in range(Ref.shape[1] - 100):
                    if Ref[i + 1][-j - 100] != 65533:
                        long_tie_length2 = j + 100
                        break

                tow_tie_sign = 1

        if tow_tie_sign == 1:
            break

    for i in range(Ref.shape[0]):
        if i > 4 and i < (Ref.shape[0] - 4):
            if Ref[i][0] == 65533 and Ref[i - 1][0] != 65533:
                for j in range(short_tie_length):
                    Ref[i][j] = Ref[i - 1][j]
                    Ref[i + 3][j] = Ref[i + 4][j]

                for j in range(short_tie_length2):
                    Ref[i][-j] = Ref[i - 1][-j]
                    Ref[i + 3][-j] = Ref[i + 4][-j]

                for j in range(long_tie_length):
                    Ref[i + 1][j] = Ref[i][j]
                    Ref[i + 2][j] = Ref[i][j]

                for j in range(long_tie_length2):
                    Ref[i + 1][-j] = Ref[i + 3][-j]
                    Ref[i + 2][-j] = Ref[i + 3][-j]

    return Ref


def viirs_i2(Ref):
    tow_tie_sign = 0
    short_tie_length = 0
    long_tie_length = 0
    short_tie_length2 = 0
    long_tie_length2 = 0
    most_tie_length = 0
    most_tie_length2 = 0
    for i in range(Ref.shape[0]):
        if i > 10:
            if Ref[i][2] == -9999:
                for j in range(Ref.shape[1] - 100):
                    if Ref[i][j + 100] != -9999:
                        short_tie_length = j + 100
                        break
                for j in range(Ref.shape[1] - 100):
                    if Ref[i][-j - 100] != -9999:
                        short_tie_length2 = j + 100
                        break
                for j in range(Ref.shape[1] - 100):
                    if Ref[i + 2][j + 100] != -9999:
                        long_tie_length = j + 100
                        break

                for j in range(Ref.shape[1] - 100):
                    if Ref[i + 2][-j - 100] != -9999:
                        long_tie_length2 = j + 100
                        break

                for j in range(Ref.shape[1] - 100):
                    if Ref[i + 4][j + 100] != -9999:
                        most_tie_length = j + 100
                        break

                for j in range(Ref.shape[1] - 100):
                    if Ref[i + 4][-j - 100] != -9999:
                        most_tie_length2 = j + 100
                        break

                tow_tie_sign = 1

        if tow_tie_sign == 1:
            break

    for i in range(Ref.shape[0]):
        if i > 4 and i < (Ref.shape[0] - 18):
            if Ref[i][0] == -9999 and Ref[i - 1][0] != -9999:
                for j in range(short_tie_length):
                    Ref[i][j] = Ref[i - 1][j]
                    Ref[i + 1][j] = Ref[i - 1][j]
                    Ref[i + 6][j] = Ref[i + 8][j]
                    Ref[i + 7][j] = Ref[i + 8][j]

                for j in range(short_tie_length2):
                    Ref[i][-j] = Ref[i - 1][-j]
                    Ref[i + 1][-j] = Ref[i - 1][-j]
                    Ref[i + 10][-j] = Ref[i + 12][-j]
                    Ref[i + 11][-j] = Ref[i + 12][-j]

                for j in range(long_tie_length):
                    Ref[i + 2][j] = Ref[i][j]
                    Ref[i + 3][j] = Ref[i][j]
                    Ref[i + 8][j] = Ref[i + 10][j]
                    Ref[i + 9][j] = Ref[i + 10][j]

                for j in range(long_tie_length2):
                    Ref[i + 2][-j] = Ref[i][-j]
                    Ref[i + 3][-j] = Ref[i][-j]
                    Ref[i + 8][-j] = Ref[i + 10][-j]
                    Ref[i + 9][-j] = Ref[i + 10][-j]

                for j in range(most_tie_length):
                    Ref[i + 4][j] = Ref[i][j]
                    Ref[i + 5][j] = Ref[i][j]
                    Ref[i + 6][j] = Ref[i + 8][j]
                    Ref[i + 7][j] = Ref[i + 9][j]

                for j in range(most_tie_length2):
                    Ref[i + 4][-j] = Ref[i][-j]
                    Ref[i + 5][-j] = Ref[i][-j]
                    Ref[i + 6][-j] = Ref[i + 8][-j]
                    Ref[i + 7][-j] = Ref[i + 9][-j]

    return Ref
