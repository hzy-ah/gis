import math
import re
import sys
import os
import basic_function as bf
import create_data as cd
import read_data as rd
import zhishu as zs
import fmask_h
import Mon_GPP_NPP_test as mgnt
import LST_h
import numpy as np
import netCDF4 as nc
from osgeo import gdal
from scipy import ndimage


def NPP_VIIRS_InsProducts(data_dir, time_type, out_dir):
    res = re.search('^[0-9]{4}-[0-9]{2}-[0-9]{2}$', time_type)
    if res is None:
        print('Please check your input date!')
        return False
    year = time_type[0:4]
    month = time_type[5:7]
    day = time_type[8:10]

    file_list = os.listdir(data_dir)
    file_list_need = list()
    for file in file_list:
        if file.find(year + month + day) != -1:
            file_list_need.append(file)

    if file_list_need == []:
        print('Can not find the idea observation at your input day!')
        return False

    file_list_gimgo = list()
    file_list_gmodo = list()

    file_list_svi01 = list()
    file_list_svi02 = list()
    file_list_svi03 = list()
    file_list_svi04 = list()
    file_list_svi05 = list()

    file_list_svm02 = list()
    file_list_svm04 = list()
    file_list_svm08 = list()
    file_list_svm09 = list()
    file_list_svm11 = list()
    file_list_svm15 = list()
    file_list_svm16 = list()
    # 根据i波段和m波段的名称进行文件匹配
    for file in file_list_need:
        if file.lower().find('svi01') != -1:
            file_list_svi01.append(file)
        if file.lower().find('svi02') != -1:
            file_list_svi02.append(file)
        if file.lower().find('svi03') != -1:
            file_list_svi03.append(file)
        if file.lower().find('svi04') != -1:
            file_list_svi04.append(file)
        if file.lower().find('svi05') != -1:
            file_list_svi05.append(file)
        if file.lower().find('svm02') != -1:
            file_list_svm02.append(file)
        if file.lower().find('svm04') != -1:
            file_list_svm04.append(file)
        if file.lower().find('svm08') != -1:
            file_list_svm08.append(file)
        if file.lower().find('svm09') != -1:
            file_list_svm09.append(file)
        if file.lower().find('svm11') != -1:
            file_list_svm11.append(file)
        if file.lower().find('svm15') != -1:
            file_list_svm15.append(file)
        if file.lower().find('svm16') != -1:
            file_list_svm16.append(file)
        if file.lower().find('gimgo') != -1:
            file_list_gimgo.append(file)
        if file.lower().find('gmodo') != -1:
            file_list_gmodo.append(file)

    Ref_all = list()
    Bt_all = list()
    Ref2_all = list()
    Bt2_all = list()
    Sa_all = list()
    So_all = list()
    A_all = list()
    Sa2_all = list()
    So2_all = list()
    A2_all = list()

    for data_num in range(len(file_list_svi01)):

        gimgo_path = data_dir + '/' + file_list_gimgo[data_num]
        gmodo_path = data_dir + '/' + file_list_gmodo[data_num]

        Ref_svi = []
        Bt_svi = []

        # 判断文件在福建省范围内是否有足够的像素点
        if bf.ifuse(gimgo_path, None, 'VIIRS') == -1:
            continue

        h1, lat1, lon1, sa_z1, so_z1, A1 = bf.geo(gimgo_path, data_type='VIIRS')
        # 设置初始值
        x_max = 0
        x_min = 12000
        y_max = 0
        y_min = 12000
        # 进行初步的裁剪，减少辐射定标的数据量
        for m in range(lat1.shape[0]):
            for n in range(lat1.shape[1]):
                if lat1[m][n] >= 20 and lat1[m][n] <= 30 and lon1[m][n] >= 114 and lon1[m][n] <= 124:
                    if m > x_max:
                        x_max = m

                    if m < x_min:
                        x_min = m

                    if n > y_max:
                        y_max = n

                    if n < y_min:
                        y_min = n

        data_range = (x_min, x_max, y_min, y_max)

        lat1 = lat1[x_min:x_max + 1, y_min:y_max + 1]
        lon1 = lon1[x_min:x_max + 1, y_min:y_max + 1]
        sa_z1 = sa_z1[x_min:x_max + 1, y_min:y_max + 1]
        so_z1 = so_z1[x_min:x_max + 1, y_min:y_max + 1]
        A1 = A1[x_min:x_max + 1, y_min:y_max + 1]

        file_path = data_dir + '/' + file_list_svi01[data_num]
        ref_svi01 = bf.band(file_path, data_type='VIIRS', band_type='ref_i', data_range=data_range)
        Ref_svi.append(ref_svi01)

        file_path = data_dir + '/' + file_list_svi02[data_num]
        ref_svi02 = bf.band(file_path, data_type='VIIRS', band_type='ref_i', data_range=data_range)
        Ref_svi.append(ref_svi02)

        file_path = data_dir + '/' + file_list_svi03[data_num]
        ref_svi03 = bf.band(file_path, data_type='VIIRS', band_type='ref_i', data_range=data_range)
        Ref_svi.append(ref_svi03)

        file_path = data_dir + '/' + file_list_svi04[data_num]
        bt_svi04 = bf.band(file_path, data_type='VIIRS', band_type='bt_i', data_range=data_range)
        Bt_svi.append(bt_svi04)

        file_path = data_dir + '/' + file_list_svi05[data_num]
        bt_svi05 = bf.band(file_path, data_type='VIIRS', band_type='bt_i', data_range=data_range)
        Bt_svi.append(bt_svi05)

        Ref_svi = np.array(Ref_svi).astype(float)
        Bt_svi = np.array(Bt_svi).astype(float)

        print(so_z1.shape)
        print(Ref_svi.shape, Bt_svi.shape)

        # 进行太阳天顶角的误差修正
        for i in range(Ref_svi.shape[1]):
            for j in range(Ref_svi.shape[2]):
                for num in range(Ref_svi.shape[0]):
                    Ref_svi[num][i][j] = Ref_svi[num][i][j] / math.cos(so_z1[i][j] * math.pi / 180)

        lat_path1 = cd.create_tif(lat1, out_dir, cd.varname(lat1))
        lon_path1 = cd.create_tif(lon1, out_dir, cd.varname(lon1))
        sa_z_path1 = cd.create_tif(sa_z1, out_dir, cd.varname(sa_z1))
        so_z_path1 = cd.create_tif(so_z1, out_dir, cd.varname(so_z1))
        A_path1 = cd.create_tif(A1, out_dir, cd.varname(A1))

        vrt_path = out_dir + '/' + 'vrt.xml'

        outputBounds = (115, 22, 122, 29)
        xRes = 0.005
        yRes = 0.005

        ref = np.zeros((Ref_svi.shape[0], 1400, 1400))
        bt = np.zeros((Bt_svi.shape[0], 1400, 1400))
        # 进行i波段的几何校正，分为反射波段和亮温波段
        for i in range(Ref_svi.shape[0]):
            ref_path1 = cd.create_tif(Ref_svi[i], out_dir, cd.varname(Ref_svi))
            cd.xml_file(ref_path1, lat_path1, lon_path1, vrt_path, Ref_svi[i].shape[0], Ref_svi[i].shape[1])
            bf.renaming(vrt_path)
            swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(Ref_svi), outputBounds, xRes, yRes)
            ref[i] = swath_data.ReadAsArray()
            ref[i] = cd.replace_fill_value(ref[i])
            swath_data = None
            os.remove(out)
        os.remove(ref_path1)

        for j in range(Bt_svi.shape[0]):
            emissive_path1 = cd.create_tif(Bt_svi[j], out_dir, cd.varname(Bt_svi))
            cd.xml_file(emissive_path1, lat_path1, lon_path1, vrt_path, Bt_svi[j].shape[0], Bt_svi[j].shape[1])
            bf.renaming(vrt_path)
            swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(Bt_svi), outputBounds, xRes, yRes)
            bt[j] = cd.replace_fill_value(swath_data.ReadAsArray())
            swath_data = None
            os.remove(out)
        os.remove(emissive_path1)

        cd.xml_file(sa_z_path1, lat_path1, lon_path1, vrt_path, sa_z1.shape[0], sa_z1.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(sa_z1), outputBounds, xRes, yRes)
        sa_z = cd.replace_fill_value(swath_data.ReadAsArray())
        swath_data = None
        os.remove(out)

        cd.xml_file(so_z_path1, lat_path1, lon_path1, vrt_path, so_z1.shape[0], so_z1.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(so_z1), outputBounds, xRes, yRes)
        so_z = cd.replace_fill_value(swath_data.ReadAsArray())
        swath_data = None
        os.remove(out)

        cd.xml_file(A_path1, lat_path1, lon_path1, vrt_path, A1.shape[0], A1.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(A1), outputBounds, xRes, yRes)
        A = cd.replace_fill_value(swath_data.ReadAsArray())
        lat_nc, lon_nc = rd.read_lat_lon(swath_data)
        swath_data = None
        os.remove(out)

        Ref_all.append(ref)
        Bt_all.append(bt)
        Sa_all.append(sa_z)
        So_all.append(so_z)
        A_all.append(A)

        os.remove(lat_path1)
        os.remove(lon_path1)
        os.remove(sa_z_path1)
        os.remove(so_z_path1)
        os.remove(A_path1)

        Ref_svm = []
        Bt_svm = []

        band_type = [2, 4, 8, 9, 11, 15, 16]

        file_path = data_dir + '/' + file_list_svm02[data_num]
        ref_svm02 = bf.band(file_path, data_type='VIIRS', band_type='ref_m')
        Ref_svm.append(ref_svm02)

        file_path = data_dir + '/' + file_list_svm04[data_num]
        ref_svm04 = bf.band(file_path, data_type='VIIRS', band_type='ref_m')
        Ref_svm.append(ref_svm04)

        file_path = data_dir + '/' + file_list_svm08[data_num]
        ref_svm08 = bf.band(file_path, data_type='VIIRS', band_type='ref_m')
        Ref_svm.append(ref_svm08)

        file_path = data_dir + '/' + file_list_svm09[data_num]
        ref_svm09 = bf.band(file_path, data_type='VIIRS', band_type='ref_m')
        Ref_svm.append(ref_svm09)

        file_path = data_dir + '/' + file_list_svm11[data_num]
        ref_svm11 = bf.band(file_path, data_type='VIIRS', band_type='ref_m')
        Ref_svm.append(ref_svm11)

        file_path = data_dir + '/' + file_list_svm15[data_num]
        bt_svm15 = bf.band(file_path, data_type='VIIRS', band_type='bt_m')
        Bt_svm.append(bt_svm15)

        file_path = data_dir + '/' + file_list_svm16[data_num]
        bt_svm16 = bf.band(file_path, data_type='VIIRS', band_type='bt_m')
        Bt_svm.append(bt_svm16)

        Ref_svm = np.array(Ref_svm).astype(float)
        Bt_svm = np.array(Bt_svm).astype(float)

        h2, lat2, lon2, sa_z2, so_z2, A2 = bf.geo(gmodo_path, data_type='VIIRS')

        # 进行太阳天顶角误差校正
        for i in range(Ref_svm.shape[1]):
            for j in range(Ref_svm.shape[2]):
                for num in range(Ref_svm.shape[0]):
                    Ref_svm[num][i][j] = Ref_svm[num][i][j] / math.cos(so_z2[i][j] * math.pi / 180)

        lat_path2 = cd.create_tif(lat2, out_dir, cd.varname(lat2))
        lon_path2 = cd.create_tif(lon2, out_dir, cd.varname(lon2))
        sa_z_path2 = cd.create_tif(sa_z2, out_dir, cd.varname(sa_z2))
        so_z_path2 = cd.create_tif(so_z2, out_dir, cd.varname(so_z2))
        A_path2 = cd.create_tif(A2, out_dir, cd.varname(A2))

        vrt_path = out_dir + '/' + 'vrt.xml'

        ref2 = np.zeros((Ref_svm.shape[0], 700, 700))
        bt2 = np.zeros((Bt_svm.shape[0], 700, 700))
        xRes = 0.01
        yRes = 0.01

        # 进行m波段的几何校正，分为反射波段和亮温波段
        for i in range(Ref_svm.shape[0]):
            ref_path2 = cd.create_tif(Ref_svm[i], out_dir, cd.varname(Ref_svm))
            cd.xml_file(ref_path2, lat_path2, lon_path2, vrt_path, Ref_svm[i].shape[0], Ref_svm[i].shape[1])
            bf.renaming(vrt_path)
            swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(Ref_svm), outputBounds, xRes, yRes)
            ref2[i] = swath_data.ReadAsArray()
            ref2[i] = cd.replace_fill_value(ref2[i])
            swath_data = None
            os.remove(out)

        for j in range(Bt_svm.shape[0]):
            emissive_path2 = cd.create_tif(Bt_svm[j], out_dir, cd.varname(Bt_svm))
            cd.xml_file(emissive_path2, lat_path2, lon_path2, vrt_path, Bt_svm[j].shape[0], Bt_svm[j].shape[1])
            bf.renaming(vrt_path)
            swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(Bt_svm), outputBounds, xRes, yRes)
            bt2[j] = cd.replace_fill_value(swath_data.ReadAsArray())
            swath_data = None
            os.remove(out)

        cd.xml_file(sa_z_path2, lat_path2, lon_path2, vrt_path, sa_z2.shape[0], sa_z2.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(sa_z2), outputBounds, xRes, yRes)
        sa_z2 = cd.replace_fill_value(swath_data.ReadAsArray())
        swath_data = None
        os.remove(out)

        cd.xml_file(so_z_path2, lat_path2, lon_path2, vrt_path, so_z2.shape[0], so_z2.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(so_z2), outputBounds, xRes, yRes)
        so_z2 = cd.replace_fill_value(swath_data.ReadAsArray())
        swath_data = None
        os.remove(out)

        cd.xml_file(A_path2, lat_path2, lon_path2, vrt_path, A2.shape[0], A2.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(A2), outputBounds, xRes, yRes)
        A2 = cd.replace_fill_value(swath_data.ReadAsArray())
        lat_nc2, lon_nc2 = rd.read_lat_lon(swath_data)
        swath_data = None
        os.remove(out)

        Ref2_all.append(ref2)
        Bt2_all.append(bt2)
        Sa2_all.append(sa_z2)
        So2_all.append(so_z2)
        A2_all.append(A2)

        os.remove(ref_path2)
        os.remove(emissive_path2)
        os.remove(lat_path2)
        os.remove(lon_path2)
        os.remove(sa_z_path2)
        os.remove(so_z_path2)
        os.remove(A_path2)

    if Ref_all == []:
        print('Can not find the idea NPP-VIIRS L1B file at your input date and time!')
        return False

    Ref_all = np.array(Ref_all).astype(float)
    Bt_all = np.array(Bt_all).astype(float)
    Sa_all = np.array(Sa_all).astype(float)
    So_all = np.array(So_all).astype(float)
    A_all = np.array(A_all).astype(float)

    for num in range(Ref_all.shape[0]):
        for num2 in range(Ref_all.shape[1]):
            for i in range(Ref_all.shape[2]):
                for j in range(Ref_all.shape[3]):
                    if ref[num2][i][j] == -999 and Ref_all[num][num2][i][j] < 1:
                        ref[num2][i][j] = Ref_all[num][num2][i][j]

    for num in range(Bt_all.shape[0]):
        for num2 in range(Bt_all.shape[1]):
            for i in range(Bt_all.shape[2]):
                for j in range(Bt_all.shape[3]):
                    if bt[num2][i][j] == -999 and 0 < Bt_all[num][num2][i][j] < 400:
                        bt[num2][i][j] = Bt_all[num][num2][i][j]

    for num in range(A_all.shape[0]):
        for i in range(A_all.shape[1]):
            for j in range(A_all.shape[2]):
                if sa_z[i][j] == -999 and 0 < Sa_all[num][i][j] < 100:
                    sa_z[i][j] = Sa_all[num][i][j]

                if so_z[i][j] == -999 and 0 < So_all[num][i][j] < 100:
                    so_z[i][j] = So_all[num][i][j]

                if A[i][j] == -999 and 0 < A_all[num][i][j] < 180:
                    A[i][j] = A_all[num][i][j]

    Ref2_all = np.array(Ref2_all).astype(float)
    Bt2_all = np.array(Bt2_all).astype(float)
    Sa2_all = np.array(Sa2_all).astype(float)
    So2_all = np.array(So2_all).astype(float)
    A2_all = np.array(A2_all).astype(float)

    for num in range(Ref2_all.shape[0]):
        for num2 in range(Ref2_all.shape[1]):
            for i in range(Ref2_all.shape[2]):
                for j in range(Ref2_all.shape[3]):
                    if ref2[num2][i][j] == -999 and 0 < Ref2_all[num][num2][i][j] < 1:
                        ref2[num2][i][j] = Ref2_all[num][num2][i][j]

    for num in range(Bt2_all.shape[0]):
        for num2 in range(Bt2_all.shape[1]):
            for i in range(Bt2_all.shape[2]):
                for j in range(Bt2_all.shape[3]):
                    if bt2[num2][i][j] == -999 and 0 < Bt2_all[num][num2][i][j] < 400:
                        bt2[num2][i][j] = Bt2_all[num][num2][i][j]

    for num in range(A2_all.shape[0]):
        for i in range(A2_all.shape[1]):
            for j in range(A2_all.shape[2]):
                if sa_z2[i][j] == -999 and 0 < Sa2_all[num][i][j] < 100:
                    sa_z2[i][j] = Sa2_all[num][i][j]

                if so_z2[i][j] == -999 and 0 < So2_all[num][i][j] < 100:
                    so_z2[i][j] = So2_all[num][i][j]

                if A2[i][j] == -999 and 0 < A2_all[num][i][j] < 180:
                    A2[i][j] = A2_all[num][i][j]

    red_cloud = ref[0]
    nir_cloud = ref[1]
    blue_cloud = cd.expand_band(ref2[0])
    green_cloud = cd.expand_band(ref2[1])
    swir1_cloud = ref[2]
    swir2_cloud = cd.expand_band(ref2[4])
    cirrus_cloud = cd.expand_band(ref2[3])
    bt_tir1_cloud = cd.expand_band(bt2[0])
    ndvi_cloud = zs.NDVI(nir_cloud, red_cloud)
    ndsi_cloud = zs.NDSI(swir1_cloud, green_cloud)
    # 进行云掩膜的计算
    cloud = fmask_h.fmask(red_cloud, nir_cloud, blue_cloud, green_cloud, swir1_cloud,
                          swir2_cloud, cirrus_cloud, bt_tir1_cloud, ndvi_cloud, ndsi_cloud)

    # 对i波段和m波段进行大气校正
    atmos_ref_all = list()
    for i in range(3):
        ref_i_whole = [1, 2, 3]
        ref_i_index = [0, 1, 2]
        table_file = './AT_LUTS/VIIRS_I' + str(ref_i_whole[i]) + '.nc'
        atmos_ref = bf.atmos(ref[int(ref_i_index[i])], sa_z, so_z, A, table_file, cloud)
        atmos_ref_all.append(atmos_ref)

    atmos_ref_all = np.array(atmos_ref_all).astype(float)

    atmos_ref2_all = list()
    for i in range(4):
        ref_i_whole = [2, 4, 8, 11]
        ref_i_index = [0, 1, 2, 4]
        table_file = './AT_LUTS/VIIRS_M' + str(ref_i_whole[i]) + '.nc'
        atmos_ref2 = bf.atmos(ref2[int(ref_i_index[i])], sa_z2, so_z2, A2, table_file, cloud[::2, ::2])
        atmos_ref2_all.append(atmos_ref2)

    atmos_ref2_all = np.array(atmos_ref2_all).astype(float)

    blue = cd.expand_band(atmos_ref2_all[0])
    red = atmos_ref_all[0]
    nir = atmos_ref_all[1]
    nir2 = atmos_ref_all[2]
    f = gdal.Open('./ndata/IGBP_lai_500.tif')
    landcover = f.ReadAsArray()

    rvi = zs.RVI(nir, red)
    ndvi = zs.NDVI(nir, red)
    savi = zs.SAVI(nir, red, 0.5)
    gvi = zs.GVI(blue, red, nir, nir2)
    evi = zs.EVI(nir, red, blue)
    lai = zs.LAI(ndvi, landcover)

    for i in range(ndvi.shape[0]):
        for j in range(ndvi.shape[1]):
            if cloud[i][j] == 1:
                rvi[i][j] = -999
                ndvi[i][j] = -999
                savi[i][j] = -999
                gvi[i][j] = -999
                evi[i][j] = -999
                lai[i][j] = -999

    # 将计算好的地表生物量瞬时产品保存为nc文件
    out_head2 = '(NPP-VIIRS)'
    cd.zhishu_nc(rvi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='RVI')
    cd.zhishu_nc(ndvi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='NDVI')
    cd.zhishu_nc(savi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='SAVI')
    cd.zhishu_nc(gvi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='GVI')
    cd.zhishu_nc(evi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='EVI')
    cd.zhishu_nc(lai, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='LAI')


def JPSS1_VIIRS_InsProducts(data_dir, time_type, out_dir):
    res = re.search('^[0-9]{4}-[0-9]{2}-[0-9]{2}$', time_type)
    if res is None:
        print('Please check your input date!')
        return False
    year = time_type[0:4]
    month = time_type[5:7]
    day = time_type[8:10]

    file_list = os.listdir(data_dir)
    file_list_need = list()
    for file in file_list:
        if file.find(year + month + day) != -1:
            file_list_need.append(file)

    if file_list_need == []:
        print('Can not find the idea observation at your input day!')
        return False

    file_list_gimgo = list()
    file_list_gmodo = list()

    file_list_svi01 = list()
    file_list_svi02 = list()
    file_list_svi03 = list()
    file_list_svi04 = list()
    file_list_svi05 = list()

    file_list_svm02 = list()
    file_list_svm04 = list()
    file_list_svm08 = list()
    file_list_svm09 = list()
    file_list_svm11 = list()
    file_list_svm15 = list()
    file_list_svm16 = list()

    for file in file_list_need:
        if file.lower().find('svi01') != -1:
            file_list_svi01.append(file)
        if file.lower().find('svi02') != -1:
            file_list_svi02.append(file)
        if file.lower().find('svi03') != -1:
            file_list_svi03.append(file)
        if file.lower().find('svi04') != -1:
            file_list_svi04.append(file)
        if file.lower().find('svi05') != -1:
            file_list_svi05.append(file)
        if file.lower().find('svm02') != -1:
            file_list_svm02.append(file)
        if file.lower().find('svm04') != -1:
            file_list_svm04.append(file)
        if file.lower().find('svm08') != -1:
            file_list_svm08.append(file)
        if file.lower().find('svm09') != -1:
            file_list_svm09.append(file)
        if file.lower().find('svm11') != -1:
            file_list_svm11.append(file)
        if file.lower().find('svm15') != -1:
            file_list_svm15.append(file)
        if file.lower().find('svm16') != -1:
            file_list_svm16.append(file)
        if file.lower().find('gimgo') != -1:
            file_list_gimgo.append(file)
        if file.lower().find('gmodo') != -1:
            file_list_gmodo.append(file)

    Ref_all = list()
    Bt_all = list()
    Ref2_all = list()
    Bt2_all = list()
    Sa_all = list()
    So_all = list()
    A_all = list()
    Sa2_all = list()
    So2_all = list()
    A2_all = list()

    for data_num in range(len(file_list_svi01)):

        gimgo_path = data_dir + '/' + file_list_gimgo[data_num]
        gmodo_path = data_dir + '/' + file_list_gmodo[data_num]

        Ref_svi = []
        Bt_svi = []

        if bf.ifuse(gimgo_path, None, 'VIIRS') == -1:
            continue

        h1, lat1, lon1, sa_z1, so_z1, A1 = bf.geo(gimgo_path, data_type='VIIRS')

        x_max = 0
        x_min = 10000
        y_max = 0
        y_min = 10000
        # 进行初步裁剪，减少辐射定标数据量
        for m in range(lat1.shape[0]):
            for n in range(lat1.shape[1]):
                if lat1[m][n] >= 20 and lat1[m][n] <= 30 and lon1[m][n] >= 114 and lon1[m][n] <= 124:
                    if m > x_max:
                        x_max = m

                    if m < x_min:
                        x_min = m

                    if n > y_max:
                        y_max = n

                    if n < y_min:
                        y_min = n

        data_range = (x_min, x_max, y_min, y_max)

        lat1 = lat1[x_min:x_max + 1, y_min:y_max + 1]
        lon1 = lon1[x_min:x_max + 1, y_min:y_max + 1]
        sa_z1 = sa_z1[x_min:x_max + 1, y_min:y_max + 1]
        so_z1 = so_z1[x_min:x_max + 1, y_min:y_max + 1]
        A1 = A1[x_min:x_max + 1, y_min:y_max + 1]

        file_path = data_dir + '/' + file_list_svi01[data_num]
        ref_svi01 = bf.band(file_path, data_type='VIIRS', band_type='ref_i', data_range=data_range)
        Ref_svi.append(ref_svi01)

        file_path = data_dir + '/' + file_list_svi02[data_num]
        ref_svi02 = bf.band(file_path, data_type='VIIRS', band_type='ref_i', data_range=data_range)
        Ref_svi.append(ref_svi02)

        file_path = data_dir + '/' + file_list_svi03[data_num]
        ref_svi03 = bf.band(file_path, data_type='VIIRS', band_type='ref_i', data_range=data_range)
        Ref_svi.append(ref_svi03)

        file_path = data_dir + '/' + file_list_svi04[data_num]
        bt_svi04 = bf.band(file_path, data_type='VIIRS', band_type='bt_i', data_range=data_range)
        Bt_svi.append(bt_svi04)

        file_path = data_dir + '/' + file_list_svi05[data_num]
        bt_svi05 = bf.band(file_path, data_type='VIIRS', band_type='bt_i', data_range=data_range)
        Bt_svi.append(bt_svi05)

        Ref_svi = np.array(Ref_svi).astype(float)
        Bt_svi = np.array(Bt_svi).astype(float)
        # 裁剪i波段缺失的条带，避免几何校正发生错误
        new_Ref_svi = []
        new_Bt_svi = []
        value_to_remove = -800
        new_so_z1, removed_rows = cd.remove_rows_with_value(so_z1, value_to_remove)
        print(removed_rows)
        new_lat1 = np.delete(lat1, removed_rows, axis=0)
        new_lon1 = np.delete(lon1, removed_rows, axis=0)
        new_sa_z1 = np.delete(sa_z1, removed_rows, axis=0)
        new_A1 = np.delete(A1, removed_rows, axis=0)
        for i in range(Ref_svi.shape[0]):
            new_ref = np.delete(Ref_svi[i], removed_rows, axis=0)
            new_Ref_svi.append(new_ref)

        for j in range(Bt_svi.shape[0]):
            new_bt = np.delete(Bt_svi[j], removed_rows, axis=0)
            new_Bt_svi.append(new_bt)

        new_Ref_svi = np.array(new_Ref_svi)
        new_Bt_svi = np.array(new_Bt_svi)

        lat_path1 = cd.create_tif(new_lat1, out_dir, cd.varname(new_lat1))
        lon_path1 = cd.create_tif(new_lon1, out_dir, cd.varname(new_lon1))
        sa_z_path1 = cd.create_tif(new_sa_z1, out_dir, cd.varname(new_sa_z1))
        so_z_path1 = cd.create_tif(new_so_z1, out_dir, cd.varname(new_so_z1))
        A_path1 = cd.create_tif(new_A1, out_dir, cd.varname(new_A1))

        vrt_path = out_dir + '/' + 'vrt.xml'

        outputBounds = (115, 22, 122, 29)
        xRes = 0.005
        yRes = 0.005

        ref = np.zeros((new_Ref_svi.shape[0], 1400, 1400))
        bt = np.zeros((new_Bt_svi.shape[0], 1400, 1400))
        # 进行i波段的几何校正
        for i in range(new_Ref_svi.shape[0]):
            ref_path1 = cd.create_tif(new_Ref_svi[i], out_dir, cd.varname(new_Ref_svi))
            cd.xml_file(ref_path1, lat_path1, lon_path1, vrt_path, new_Ref_svi[i].shape[0], new_Ref_svi[i].shape[1])
            bf.renaming(vrt_path)
            swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(new_Ref_svi), outputBounds, xRes, yRes)
            ref[i] = swath_data.ReadAsArray()
            ref[i] = cd.replace_fill_value(ref[i])
            swath_data = None
            os.remove(out)
        os.remove(ref_path1)

        for j in range(new_Bt_svi.shape[0]):
            emissive_path1 = cd.create_tif(new_Bt_svi[j], out_dir, cd.varname(new_Bt_svi))
            cd.xml_file(emissive_path1, lat_path1, lon_path1, vrt_path, new_Bt_svi[j].shape[0], new_Bt_svi[j].shape[1])
            bf.renaming(vrt_path)
            swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(new_Bt_svi), outputBounds, xRes, yRes)
            bt[j] = cd.replace_fill_value(swath_data.ReadAsArray())
            swath_data = None
            os.remove(out)
        os.remove(emissive_path1)

        cd.xml_file(sa_z_path1, lat_path1, lon_path1, vrt_path, sa_z1.shape[0], sa_z1.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(sa_z1), outputBounds, xRes, yRes)
        sa_z = cd.replace_fill_value(swath_data.ReadAsArray())
        swath_data = None
        os.remove(out)

        cd.xml_file(so_z_path1, lat_path1, lon_path1, vrt_path, so_z1.shape[0], so_z1.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(so_z1), outputBounds, xRes, yRes)
        so_z = cd.replace_fill_value(swath_data.ReadAsArray())
        swath_data = None
        os.remove(out)

        cd.xml_file(A_path1, lat_path1, lon_path1, vrt_path, A1.shape[0], A1.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(A1), outputBounds, xRes, yRes)
        A = cd.replace_fill_value(swath_data.ReadAsArray())
        lat_nc, lon_nc = rd.read_lat_lon(swath_data)
        swath_data = None
        os.remove(out)

        Ref_all.append(ref)
        Bt_all.append(bt)
        Sa_all.append(sa_z)
        So_all.append(so_z)
        A_all.append(A)

        os.remove(lat_path1)
        os.remove(lon_path1)
        os.remove(sa_z_path1)
        os.remove(so_z_path1)
        os.remove(A_path1)

        Ref_svm = []
        Bt_svm = []

        band_type = [2, 4, 8, 9, 11, 15, 16]

        file_path = data_dir + '/' + file_list_svm02[data_num]
        ref_svm02 = bf.band(file_path, data_type='VIIRS', band_type='ref_m')
        Ref_svm.append(ref_svm02)

        file_path = data_dir + '/' + file_list_svm04[data_num]
        ref_svm04 = bf.band(file_path, data_type='VIIRS', band_type='ref_m')
        Ref_svm.append(ref_svm04)

        file_path = data_dir + '/' + file_list_svm08[data_num]
        ref_svm08 = bf.band(file_path, data_type='VIIRS', band_type='ref_m')
        Ref_svm.append(ref_svm08)

        file_path = data_dir + '/' + file_list_svm09[data_num]
        ref_svm09 = bf.band(file_path, data_type='VIIRS', band_type='ref_m')
        Ref_svm.append(ref_svm09)

        file_path = data_dir + '/' + file_list_svm11[data_num]
        ref_svm11 = bf.band(file_path, data_type='VIIRS', band_type='ref_m')
        Ref_svm.append(ref_svm11)

        file_path = data_dir + '/' + file_list_svm15[data_num]

        bt_svm15 = bf.band(file_path, data_type='VIIRS', band_type='bt_m')
        Bt_svm.append(bt_svm15)

        file_path = data_dir + '/' + file_list_svm16[data_num]

        bt_svm16 = bf.band(file_path, data_type='VIIRS', band_type='bt_m')
        Bt_svm.append(bt_svm16)

        Ref_svm = np.array(Ref_svm).astype(float)
        Bt_svm = np.array(Bt_svm).astype(float)

        print(Ref_svm, Bt_svm)

        h2, lat2, lon2, sa_z2, so_z2, A2 = bf.geo(gmodo_path, data_type='VIIRS')
        # 裁剪m波段缺失的条带，避免几何校正发生错误
        new_Ref_svm = []
        new_Bt_svm = []
        value_to_remove = -800
        new_so_z2, removed_rows2 = cd.remove_rows_with_value(so_z2, value_to_remove)
        new_lat2 = np.delete(lat2, removed_rows2, axis=0)
        new_lon2 = np.delete(lon2, removed_rows2, axis=0)
        new_sa_z2 = np.delete(sa_z2, removed_rows2, axis=0)
        new_A2 = np.delete(A2, removed_rows2, axis=0)
        for i in range(Ref_svm.shape[0]):
            new_ref = np.delete(Ref_svm[i], removed_rows2, axis=0)
            new_Ref_svm.append(new_ref)

        for j in range(Bt_svm.shape[0]):
            new_bt = np.delete(Bt_svm[j], removed_rows2, axis=0)
            new_Bt_svm.append(new_bt)

        new_Ref_svm = np.array(new_Ref_svm)
        new_Bt_svm = np.array(new_Bt_svm)

        lat_path2 = cd.create_tif(new_lat2, out_dir, cd.varname(new_lat2))
        lon_path2 = cd.create_tif(new_lon2, out_dir, cd.varname(new_lon2))
        sa_z_path2 = cd.create_tif(new_sa_z2, out_dir, cd.varname(new_sa_z2))
        so_z_path2 = cd.create_tif(new_so_z2, out_dir, cd.varname(new_so_z2))
        A_path2 = cd.create_tif(new_A2, out_dir, cd.varname(new_A2))

        vrt_path = out_dir + '/' + 'vrt.xml'

        ref2 = np.zeros((Ref_svm.shape[0], 700, 700))
        bt2 = np.zeros((Bt_svm.shape[0], 700, 700))
        xRes = 0.01
        yRes = 0.01
        # 进行m波段的几何校正
        for i in range(new_Ref_svm.shape[0]):
            ref_path2 = cd.create_tif(new_Ref_svm[i], out_dir, cd.varname(new_Ref_svm))
            cd.xml_file(ref_path2, lat_path2, lon_path2, vrt_path, new_Ref_svm[i].shape[0], new_Ref_svm[i].shape[1])
            bf.renaming(vrt_path)
            swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(new_Ref_svm), outputBounds, xRes, yRes)
            ref2[i] = swath_data.ReadAsArray()
            ref2[i] = cd.replace_fill_value(ref2[i])
            swath_data = None
            os.remove(out)

        for j in range(new_Bt_svm.shape[0]):
            emissive_path2 = cd.create_tif(new_Bt_svm[j], out_dir, cd.varname(new_Bt_svm))
            cd.xml_file(emissive_path2, lat_path2, lon_path2, vrt_path, new_Bt_svm[j].shape[0], new_Bt_svm[j].shape[1])
            bf.renaming(vrt_path)
            swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(new_Bt_svm), outputBounds, xRes, yRes)
            bt2[j] = cd.replace_fill_value(swath_data.ReadAsArray())
            swath_data = None
            os.remove(out)

        cd.xml_file(sa_z_path2, lat_path2, lon_path2, vrt_path, sa_z2.shape[0], sa_z2.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(sa_z2), outputBounds, xRes, yRes)
        sa_z2 = cd.replace_fill_value(swath_data.ReadAsArray())
        swath_data = None
        os.remove(out)

        cd.xml_file(so_z_path2, lat_path2, lon_path2, vrt_path, so_z2.shape[0], so_z2.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(so_z2), outputBounds, xRes, yRes)
        so_z2 = cd.replace_fill_value(swath_data.ReadAsArray())
        swath_data = None
        os.remove(out)

        cd.xml_file(A_path2, lat_path2, lon_path2, vrt_path, A2.shape[0], A2.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(A2), outputBounds, xRes, yRes)
        A2 = cd.replace_fill_value(swath_data.ReadAsArray())
        lat_nc2, lon_nc2 = rd.read_lat_lon(swath_data)
        swath_data = None
        os.remove(out)

        Ref2_all.append(ref2)
        Bt2_all.append(bt2)
        Sa2_all.append(sa_z2)
        So2_all.append(so_z2)
        A2_all.append(A2)

        os.remove(ref_path2)
        os.remove(emissive_path2)
        os.remove(lat_path2)
        os.remove(lon_path2)
        os.remove(sa_z_path2)
        os.remove(so_z_path2)
        os.remove(A_path2)

    if Ref_all == []:
        print('Can not find the idea JPSS1-VIIRS L1B file at your input date and time!')
        return False

    Ref_all = np.array(Ref_all).astype(float)
    Bt_all = np.array(Bt_all).astype(float)
    Sa_all = np.array(Sa_all).astype(float)
    So_all = np.array(So_all).astype(float)
    A_all = np.array(A_all).astype(float)

    for num in range(Ref_all.shape[0]):
        for num2 in range(Ref_all.shape[1]):
            for i in range(Ref_all.shape[2]):
                for j in range(Ref_all.shape[3]):
                    if ref[num2][i][j] == -999 and 0 < Ref_all[num][num2][i][j] < 1:
                        ref[num2][i][j] = Ref_all[num][num2][i][j]

    for num in range(Bt_all.shape[0]):
        for num2 in range(Bt_all.shape[1]):
            for i in range(Bt_all.shape[2]):
                for j in range(Bt_all.shape[3]):
                    if bt[num2][i][j] == -999 and 0 < Bt_all[num][num2][i][j] < 400:
                        bt[num2][i][j] = Bt_all[num][num2][i][j]

    for num in range(A_all.shape[0]):
        for i in range(A_all.shape[1]):
            for j in range(A_all.shape[2]):
                if sa_z[i][j] == -999 and 0 < Sa_all[num][i][j] < 100:
                    sa_z[i][j] = Sa_all[num][i][j]

                if so_z[i][j] == -999 and 0 < So_all[num][i][j] < 100:
                    so_z[i][j] = So_all[num][i][j]

                if A[i][j] == -999 and 0 < A_all[num][i][j] < 180:
                    A[i][j] = A_all[num][i][j]

    Ref2_all = np.array(Ref2_all).astype(float)
    Bt2_all = np.array(Bt2_all).astype(float)
    Sa2_all = np.array(Sa2_all).astype(float)
    So2_all = np.array(So2_all).astype(float)
    A2_all = np.array(A2_all).astype(float)

    for num in range(Ref2_all.shape[0]):
        for num2 in range(Ref2_all.shape[1]):
            for i in range(Ref2_all.shape[2]):
                for j in range(Ref2_all.shape[3]):
                    if ref2[num2][i][j] == -999 and 0 < Ref2_all[num][num2][i][j] < 1:
                        ref2[num2][i][j] = Ref2_all[num][num2][i][j]

    for num in range(Bt2_all.shape[0]):
        for num2 in range(Bt2_all.shape[1]):
            for i in range(Bt2_all.shape[2]):
                for j in range(Bt2_all.shape[3]):
                    if bt2[num2][i][j] == -999 and 0 < Bt2_all[num][num2][i][j] < 400:
                        bt2[num2][i][j] = Bt2_all[num][num2][i][j]

    for num in range(A2_all.shape[0]):
        for i in range(A2_all.shape[1]):
            for j in range(A2_all.shape[2]):
                if sa_z2[i][j] == -999 and 0 < Sa2_all[num][i][j] < 100:
                    sa_z2[i][j] = Sa2_all[num][i][j]

                if so_z2[i][j] == -999 and 0 < So2_all[num][i][j] < 100:
                    so_z2[i][j] = So2_all[num][i][j]

                if A2[i][j] == -999 and 0 < A2_all[num][i][j] < 180:
                    A2[i][j] = A2_all[num][i][j]

    red_cloud = ref[0]
    nir_cloud = ref[1]
    blue_cloud = cd.expand_band(ref2[0])
    green_cloud = cd.expand_band(ref2[1])
    swir1_cloud = ref[2]
    swir2_cloud = cd.expand_band(ref2[4])
    cirrus_cloud = cd.expand_band(ref2[3])
    bt_tir1_cloud = cd.expand_band(bt2[0])
    ndvi_cloud = zs.NDVI(nir_cloud, red_cloud)
    ndsi_cloud = zs.NDSI(swir1_cloud, green_cloud)

    cloud = fmask_h.fmask2(nir_cloud, red_cloud, swir1_cloud, ndvi_cloud)

    # 对i波段和m波段进行大气校正
    atmos_ref_all = list()
    for i in range(3):
        ref_i_whole = [1, 2, 3]
        ref_i_index = [0, 1, 2]
        table_file = './AT_LUTS/VIIRS_I' + str(ref_i_whole[i]) + '.nc'
        atmos_ref = bf.atmos(ref[int(ref_i_index[i])], sa_z, so_z, A, table_file, cloud)
        atmos_ref_all.append(atmos_ref)

    atmos_ref_all = np.array(atmos_ref_all).astype(float)

    atmos_ref2_all = list()
    for i in range(4):
        ref_i_whole = [2, 4, 8, 11]
        ref_i_index = [0, 1, 2, 4]
        table_file = './AT_LUTS/VIIRS_M' + str(ref_i_whole[i]) + '.nc'
        atmos_ref2 = bf.atmos(ref2[int(ref_i_index[i])], sa_z2, so_z2, A2, table_file, cloud[::2, ::2])
        atmos_ref2_all.append(atmos_ref2)

    atmos_ref2_all = np.array(atmos_ref2_all).astype(float)

    blue = cd.expand_band(atmos_ref2_all[0])
    red = atmos_ref_all[0]
    nir = atmos_ref_all[1]
    nir2 = atmos_ref_all[2]
    f = gdal.Open('./ndata/IGBP_lai_500.tif')
    landcover = f.ReadAsArray()

    rvi = zs.RVI(nir, red)
    ndvi = zs.NDVI(nir, red)
    savi = zs.SAVI(nir, red, 0.5)
    gvi = zs.GVI(blue, red, nir, nir2)
    evi = zs.EVI(nir, red, blue)
    lai = zs.LAI(ndvi, landcover)

    for i in range(ndvi.shape[0]):
        for j in range(ndvi.shape[1]):
            if cloud[i][j] == 1:
                rvi[i][j] = -999
                ndvi[i][j] = -999
                savi[i][j] = -999
                gvi[i][j] = -999
                evi[i][j] = -999
                lai[i][j] = -999

    # 将地表生物量瞬时产品保存为nc文件
    out_head2 = '(JPSS1-VIIRS)'
    cd.zhishu_nc(rvi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='RVI')
    cd.zhishu_nc(ndvi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='NDVI')
    cd.zhishu_nc(savi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='SAVI')
    cd.zhishu_nc(gvi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='GVI')
    cd.zhishu_nc(evi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='EVI')
    cd.zhishu_nc(lai, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='LAI')


def FY3C_VIRR_InsProducts(data_dir, time_type, out_dir):
    res = re.search('^[0-9]{4}-[0-9]{2}-[0-9]{2}$', time_type)
    if res is None:
        print('Please check your input date!')
        return False

    year = time_type[0:4]
    month = time_type[5:7]
    day = time_type[8:10]

    file_list = os.listdir(data_dir)
    file_list_need = list()
    for file in file_list:
        if file.find(year + '_' + month + '_' + day) != -1:
            file_list_need.append(file)

    if file_list_need == []:
        print('Can not find the idea observation at your input day!')
        return False
    file_list_geo = list()
    file_list_l1b = list()

    for file in file_list_need:
        if file.find('GEO') != -1:
            file_list_geo.append(file)
        if file.find('L1B') != -1:
            file_list_l1b.append(file)

    Ref_all = list()
    Bt_all = list()
    Sa_all = list()
    So_all = list()
    A_all = list()
    for i in range(len(file_list_geo)):
        band_path = data_dir + '/' + file_list_l1b[i]
        geo_path = data_dir + '/' + file_list_geo[i]

        if bf.ifuse(geo_path, band_path, 'FY3C') == -1:
            continue
        # 读取波段、经纬度等数据
        Ref, Emissive = bf.band(band_path, data_type='FY3C')
        h, lat, lon, sa_z, so_z, A = bf.geo(geo_path, data_type='FY3C')

        lat_path = cd.create_tif(lat, out_dir, cd.varname(lat))
        lon_path = cd.create_tif(lon, out_dir, cd.varname(lon))
        sa_z_path = cd.create_tif(sa_z, out_dir, cd.varname(sa_z))
        so_z_path = cd.create_tif(so_z, out_dir, cd.varname(so_z))
        A_path = cd.create_tif(A, out_dir, cd.varname(A))

        vrt_path = out_dir + '/' + 'vrt.xml'

        outputBounds = (115, 22, 122, 29)
        xRes = 0.01
        yRes = 0.01
        band_size = int(7 / xRes)
        ref = np.zeros((Ref.shape[0], band_size, band_size))
        emissive = np.zeros((Emissive.shape[0], band_size, band_size))
        # 进行几何校正
        for i in range(Ref.shape[0]):
            ref_path = cd.create_tif(Ref[i], out_dir, cd.varname(Ref))
            cd.xml_file(ref_path, lat_path, lon_path, vrt_path, Ref[i].shape[0], Ref[i].shape[1])
            bf.renaming(vrt_path)
            swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(Ref), outputBounds, xRes, yRes)
            ref[i] = cd.replace_fill_value(swath_data.ReadAsArray())
            swath_data = None
            os.remove(out)

        for j in range(Emissive.shape[0]):
            emissive_path = cd.create_tif(Emissive[j], out_dir, cd.varname(Emissive))
            cd.xml_file(emissive_path, lat_path, lon_path, vrt_path, Emissive[j].shape[0], Emissive[j].shape[1])
            bf.renaming(vrt_path)
            swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(Emissive), outputBounds, xRes, yRes)
            emissive[j] = cd.replace_fill_value(swath_data.ReadAsArray())
            swath_data = None
            os.remove(out)

        cd.xml_file(sa_z_path, lat_path, lon_path, vrt_path, sa_z.shape[0], sa_z.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(sa_z), outputBounds, xRes, yRes)
        sa_z = cd.replace_fill_value(swath_data.ReadAsArray())
        swath_data = None
        os.remove(out)

        cd.xml_file(so_z_path, lat_path, lon_path, vrt_path, so_z.shape[0], so_z.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(so_z), outputBounds, xRes, yRes)
        so_z = cd.replace_fill_value(swath_data.ReadAsArray())
        swath_data = None
        os.remove(out)

        cd.xml_file(A_path, lat_path, lon_path, vrt_path, A.shape[0], A.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(A), outputBounds, xRes, yRes)
        A = cd.replace_fill_value(swath_data.ReadAsArray())
        lat_nc, lon_nc = rd.read_lat_lon(swath_data)
        swath_data = None
        os.remove(out)

        os.remove(ref_path)
        os.remove(emissive_path)
        os.remove(lat_path)
        os.remove(lon_path)
        os.remove(sa_z_path)
        os.remove(so_z_path)
        os.remove(A_path)

        Ref_all.append(ref)
        Bt_all.append(emissive)
        Sa_all.append(sa_z)
        So_all.append(so_z)
        A_all.append(A)

    if Ref_all == []:
        print('Can not find the idea FY3C-VIRR L1B file at your input date and time!')
        return False

    Ref_all = np.array(Ref_all).astype(float)
    Bt_all = np.array(Bt_all).astype(float)
    Sa_all = np.array(Sa_all).astype(float)
    So_all = np.array(So_all).astype(float)
    A_all = np.array(A_all).astype(float)

    for num in range(Ref_all.shape[0]):
        for num2 in range(Ref_all.shape[1]):
            for i in range(Ref_all.shape[2]):
                for j in range(Ref_all.shape[3]):
                    if ref[num2][i][j] == -999 and 0 < Ref_all[num][num2][i][j] < 1:
                        ref[num2][i][j] = Ref_all[num][num2][i][j]

    for num in range(Bt_all.shape[0]):
        for num2 in range(Bt_all.shape[1]):
            for i in range(Bt_all.shape[2]):
                for j in range(Bt_all.shape[3]):
                    if emissive[num2][i][j] == -999 and 0 < Bt_all[num][num2][i][j] < 400:
                        emissive[num2][i][j] = Bt_all[num][num2][i][j]

    for num in range(A_all.shape[0]):
        for i in range(A_all.shape[1]):
            for j in range(A_all.shape[2]):
                if sa_z[i][j] == -999 and 0 < Sa_all[num][i][j] < 100:
                    sa_z[i][j] = Sa_all[num][i][j]

                if so_z[i][j] == -999 and 0 < So_all[num][i][j] < 100:
                    so_z[i][j] = So_all[num][i][j]

                if A[i][j] == -999 and 0 < A_all[num][i][j] < 180:
                    A[i][j] = A_all[num][i][j]

    red_cloud = ref[0]
    nir_cloud = ref[1]
    blue_cloud = ref[3]
    green_cloud = ref[5]
    swir1_cloud = ref[2]
    swir2_cloud = emissive[0]
    cirrus_cloud = ref[6]
    bt_tir1_cloud = emissive[1]
    ndvi_cloud = zs.NDVI(nir_cloud, red_cloud)
    ndsi_cloud = zs.NDSI(swir1_cloud, green_cloud)
    # 计算云掩膜
    cloud = fmask_h.fmask(red_cloud, nir_cloud, blue_cloud, green_cloud, swir1_cloud,
                          swir2_cloud, cirrus_cloud, bt_tir1_cloud, ndvi_cloud, ndsi_cloud)

    # 进行大气校正
    atmos_ref_all = list()
    for i in range(5):
        ref_whole = [1, 2, 6, 7, 9]
        ref_index = [0, 1, 2, 3, 5]
        table_file = './AT_LUTS/VIRR_B' + str(ref_whole[i]) + '.nc'
        atmos_ref = bf.atmos(ref[int(ref_index[i])], sa_z, so_z, A, table_file, cloud)
        atmos_ref_all.append(atmos_ref)

    atmos_ref_all = np.array(atmos_ref_all).astype(float)

    blue = atmos_ref_all[3]
    red = atmos_ref_all[0]
    nir = atmos_ref_all[2]
    nir2 = atmos_ref_all[1]
    f = gdal.Open('./ndata/IGBP_lai_1000.tif')
    landcover = f.ReadAsArray()

    rvi = zs.RVI(nir, red)
    ndvi = zs.NDVI(nir, red)
    savi = zs.SAVI(nir, red, 0.5)
    gvi = zs.GVI(blue, red, nir2, nir)
    evi = zs.EVI(nir, red, blue)
    lai = zs.LAI(ndvi, landcover)

    for i in range(rvi.shape[0]):
        for j in range(rvi.shape[1]):
            if cloud[i][j] == 1:
                rvi[i][j] = -999
                ndvi[i][j] = -999
                savi[i][j] = -999
                gvi[i][j] = -999
                evi[i][j] = -999
                lai[i][j] = -999
    # 将地表生物量瞬时产品保存为nc文件
    out_head2 = '(FY3C-VIRR)'
    cd.zhishu_nc(rvi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='RVI')
    cd.zhishu_nc(ndvi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='NDVI')
    cd.zhishu_nc(savi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='SAVI')
    cd.zhishu_nc(gvi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='GVI')
    cd.zhishu_nc(evi, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='EVI')
    cd.zhishu_nc(lai, lat_nc, lon_nc, out_dir, year, month, day, out_head2, name='LAI')


def FY3D_MERSI_InsProducts(data_dir, time_type, out_dir):
    res = re.search('^[0-9]{4}-[0-9]{2}-[0-9]{2}$', time_type)
    if res is None:
        print('Please check your input date!')
        return False

    year = time_type[0:4]
    month = time_type[5:7]
    day = time_type[8:10]

    file_list = os.listdir(data_dir)
    file_list_need = list()
    for file in file_list:
        if file.find(year + '_' + month + '_' + day) != -1:
            file_list_need.append(file)

    if file_list_need == []:
        print('Can not find the idea observation at your input day!')
        return False

    file_list_geo_1000 = list()
    file_list_geo_250 = list()
    file_list_l1b_1000 = list()
    file_list_l1b_250 = list()

    for file in file_list_need:
        if file.find('GEO1K') != -1:
            file_list_geo_1000.append(file)
        if file.find('GEOQK') != -1:
            file_list_geo_250.append(file)

        if file.find('1000M') != -1:
            file_list_l1b_1000.append(file)

        if file.find('0250M') != -1:
            file_list_l1b_250.append(file)

    Ref_1000_all = list()
    Ref_500_all = list()
    Bt_500_all = list()
    Sa_1000_all = list()
    So_1000_all = list()
    A_1000_all = list()

    for i in range(len(file_list_geo_1000)):
        band_path_1000 = data_dir + '/' + file_list_l1b_1000[i]
        band_path_250 = data_dir + '/' + file_list_l1b_250[i]
        geo_path_1000 = data_dir + '/' + file_list_geo_1000[i]
        geo_path_250 = data_dir + '/' + file_list_geo_250[i]

        if bf.ifuse(geo_path_1000, band_path_1000, 'FY3D') == -1:
            continue
        # 读取波段、经纬度等数据
        Ref_1000 = bf.band(band_path_1000, data_type='FY3D_1000')
        h, lat_1000, lon_1000, sa_z_1000, so_z_1000, A_1000 = bf.geo(geo_path_1000, data_type='FY3D_1000')

        lat_250, lon_250, data_range_250 = bf.geo(geo_path_250, data_type='FY3D_250')

        Ref_250, Emissive_250, data_range = bf.band(band_path_250, data_type='FY3D_250', data_range=data_range_250)

        x_min = data_range[0]
        x_max = data_range[1]
        y_min = data_range[2]
        y_max = data_range[3]

        for num in range(Ref_1000.shape[0]):
            for i in range(Ref_1000.shape[1]):
                for j in range(Ref_1000.shape[2]):
                    Ref_1000[num][i][j] = Ref_1000[num][i][j] / math.cos(so_z_1000[i][j] * math.pi / 180)

        for num in range(Ref_250.shape[0]):
            for i in range(Ref_250.shape[1]):
                for j in range(Ref_250.shape[2]):
                    Ref_250[num][i][j] = Ref_250[num][i][j] / \
                                         math.cos(so_z_1000[int((i + x_min) / 4)][int((j + y_min) / 4)] * math.pi / 180)

        lat_path_250 = cd.create_tif(lat_250, out_dir, cd.varname(lat_250))
        lon_path_250 = cd.create_tif(lon_250, out_dir, cd.varname(lon_250))

        vrt_path = out_dir + '/' + 'vrt.xml'

        outputBounds = (115, 22, 122, 29)
        xRes = 0.01
        yRes = 0.01

        ref_1000 = np.zeros((Ref_1000.shape[0], 700, 700))

        x_max = 0
        x_min = 10000
        y_max = 0
        y_min = 10000
        # 进行初步裁剪，减少辐射定标工作量
        for i in range(lat_1000.shape[0]):
            for j in range(lat_1000.shape[1]):
                if lat_1000[i][j] >= 20 and lat_1000[i][j] <= 30 and lon_1000[i][j] >= 114 and lon_1000[i][j] <= 124:
                    if i > x_max:
                        x_max = i

                    if i < x_min:
                        x_min = i

                    if j > y_max:
                        y_max = j

                    if j < y_min:
                        y_min = j

        Ref_1000 = Ref_1000[:, x_min:x_max + 1, y_min:y_max + 1]
        lat_1000 = lat_1000[x_min:x_max + 1, y_min:y_max + 1]
        lon_1000 = lon_1000[x_min:x_max + 1, y_min:y_max + 1]
        sa_z_1000 = sa_z_1000[x_min:x_max + 1, y_min:y_max + 1]
        so_z_1000 = so_z_1000[x_min:x_max + 1, y_min:y_max + 1]
        A_1000 = A_1000[x_min:x_max + 1, y_min:y_max + 1]

        lat_path_1000 = cd.create_tif(lat_1000, out_dir, cd.varname(lat_1000))
        lon_path_1000 = cd.create_tif(lon_1000, out_dir, cd.varname(lon_1000))
        sa_z_path_1000 = cd.create_tif(sa_z_1000, out_dir, cd.varname(sa_z_1000))
        so_z_path_1000 = cd.create_tif(so_z_1000, out_dir, cd.varname(so_z_1000))
        A_path_1000 = cd.create_tif(A_1000, out_dir, cd.varname(A_1000))
        # 进行几何校正
        for i in range(Ref_1000.shape[0]):
            ref_path_1000 = cd.create_tif(Ref_1000[i], out_dir, cd.varname(Ref_1000))
            cd.xml_file(ref_path_1000, lat_path_1000, lon_path_1000, vrt_path, Ref_1000[i].shape[0],
                        Ref_1000[i].shape[1])
            bf.renaming(vrt_path)
            swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(Ref_1000), outputBounds, xRes, yRes)
            ref_1000[i] = cd.replace_fill_value(swath_data.ReadAsArray())
            swath_data = None
            os.remove(out)

        xRes2 = 0.005
        yRes2 = 0.005

        ref_500 = np.zeros((Ref_250.shape[0], 1400, 1400))
        emissive_500 = np.zeros((Emissive_250.shape[0], 1400, 1400))
        # 进行几何校正
        for i in range(Ref_250.shape[0]):
            ref_path_250 = cd.create_tif(Ref_250[i], out_dir, cd.varname(Ref_250))
            cd.xml_file(ref_path_250, lat_path_250, lon_path_250, vrt_path, Ref_250[i].shape[0], Ref_250[i].shape[1])
            bf.renaming(vrt_path)
            swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(Ref_250), outputBounds, xRes2, yRes2)
            ref_500[i] = cd.replace_fill_value(swath_data.ReadAsArray())
            lat_nc_500, lon_nc_500 = rd.read_lat_lon(swath_data)
            swath_data = None
            os.remove(out)
        # 进行几何校正
        for j in range(Emissive_250.shape[0]):
            emissive_path_250 = cd.create_tif(Emissive_250[j], out_dir, cd.varname(Emissive_250))
            cd.xml_file(emissive_path_250, lat_path_250, lon_path_250, vrt_path, Emissive_250[j].shape[0],
                        Emissive_250[j].shape[1])
            bf.renaming(vrt_path)
            swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(Emissive_250), outputBounds, xRes2, yRes2)
            emissive_500[j] = cd.replace_fill_value(swath_data.ReadAsArray())
            swath_data = None
            os.remove(out)

        cd.xml_file(sa_z_path_1000, lat_path_1000, lon_path_1000, vrt_path, sa_z_1000.shape[0], sa_z_1000.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(sa_z_1000), outputBounds, xRes, yRes)
        sa_z_1000 = cd.replace_fill_value(swath_data.ReadAsArray())
        swath_data = None
        os.remove(out)

        cd.xml_file(so_z_path_1000, lat_path_1000, lon_path_1000, vrt_path, so_z_1000.shape[0], so_z_1000.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(so_z_1000), outputBounds, xRes, yRes)
        so_z_1000 = cd.replace_fill_value(swath_data.ReadAsArray())
        swath_data = None
        os.remove(out)

        cd.xml_file(A_path_1000, lat_path_1000, lon_path_1000, vrt_path, A_1000.shape[0], A_1000.shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(A_1000), outputBounds, xRes, yRes)
        A_1000 = cd.replace_fill_value(swath_data.ReadAsArray())
        lat_nc_1000, lon_nc_1000 = rd.read_lat_lon(swath_data)
        swath_data = None
        os.remove(out)

        os.remove(ref_path_1000)
        os.remove(ref_path_250)
        os.remove(emissive_path_250)
        os.remove(lat_path_1000)
        os.remove(lon_path_1000)
        os.remove(lat_path_250)
        os.remove(lon_path_250)
        os.remove(sa_z_path_1000)
        os.remove(so_z_path_1000)
        os.remove(A_path_1000)

        Ref_1000_all.append(ref_1000)
        Ref_500_all.append(ref_500)
        Bt_500_all.append(emissive_500)
        Sa_1000_all.append(sa_z_1000)
        So_1000_all.append(so_z_1000)
        A_1000_all.append(A_1000)

    if Ref_1000_all == []:
        print('Can not find the idea FY3D-MERSI L1B file at your input date and time!')
        return False
    Ref_1000_all = np.array(Ref_1000_all).astype(float)
    Ref_500_all = np.array(Ref_500_all).astype(float)
    Bt_500_all = np.array(Bt_500_all).astype(float)
    Sa_1000_all = np.array(Sa_1000_all).astype(float)
    So_1000_all = np.array(So_1000_all).astype(float)
    A_1000_all = np.array(A_1000_all).astype(float)

    for num in range(Ref_1000_all.shape[0]):
        for num2 in range(Ref_1000_all.shape[1]):
            for i in range(Ref_1000_all.shape[2]):
                for j in range(Ref_1000_all.shape[3]):
                    if ref_1000[num2][i][j] == -999 and 0 < Ref_1000_all[num][num2][i][j] < 1:
                        ref_1000[num2][i][j] = Ref_1000_all[num][num2][i][j]

    for num in range(Ref_500_all.shape[0]):
        for num2 in range(Ref_500_all.shape[1]):
            for i in range(Ref_500_all.shape[2]):
                for j in range(Ref_500_all.shape[3]):
                    if ref_500[num2][i][j] == -999 and 0 < Ref_500_all[num][num2][i][j] < 1:
                        ref_500[num2][i][j] = Ref_500_all[num][num2][i][j]

    for num in range(Bt_500_all.shape[0]):
        for num2 in range(Bt_500_all.shape[1]):
            for i in range(Bt_500_all.shape[2]):
                for j in range(Bt_500_all.shape[3]):
                    if emissive_500[num2][i][j] == -999 and 0 < Bt_500_all[num][num2][i][j] < 400:
                        emissive_500[num2][i][j] = Bt_500_all[num][num2][i][j]

    for num in range(A_1000_all.shape[0]):
        for i in range(A_1000_all.shape[1]):
            for j in range(A_1000_all.shape[2]):
                if sa_z_1000[i][j] == -999 and 0 < Sa_1000_all[num][i][j] < 100:
                    sa_z_1000[i][j] = Sa_1000_all[num][i][j]

                if so_z_1000[i][j] == -999 and 0 < So_1000_all[num][i][j] < 100:
                    so_z_1000[i][j] = So_1000_all[num][i][j]

                if A_1000[i][j] == -999 and 0 < A_1000_all[num][i][j] < 180:
                    A_1000[i][j] = A_1000_all[num][i][j]

    sa_z_500 = np.zeros((1400, 1400))
    so_z_500 = np.zeros((1400, 1400))
    A_500 = np.zeros((1400, 1400))

    for i in range(1400):
        for j in range(1400):
            sa_z_500[i][j] = sa_z_1000[int(i / 2)][int(j / 2)]
            so_z_500[i][j] = so_z_1000[int(i / 2)][int(j / 2)]
            A_500[i][j] = A_1000[int(i / 2)][int(j / 2)]

    red_cloud = ref_500[2]
    nir_cloud = ref_500[3]
    swir1_cloud = cd.expand_band(ref_1000[1])
    ndvi_cloud = zs.NDVI(nir_cloud, red_cloud)

    cloud = fmask_h.fmask2(nir_cloud, red_cloud, swir1_cloud, ndvi_cloud)

    # 进行大气校正
    atmos_ref_all_500 = list()
    atmos_ref_all_1000 = list()
    for i in range(6):
        ref_whole = [1, 2, 3, 4, 6, 7]
        table_file = './AT_LUTS/MERSI_B' + str(ref_whole[i]) + '.nc'
        if ref_whole[i] < 5:
            atmos_ref = bf.atmos(ref_500[i], sa_z_500, so_z_500, A_500, table_file, cloud)
            atmos_ref_all_500.append(atmos_ref)
        if ref_whole[i] > 5:
            atmos_ref2 = bf.atmos(ref_1000[i - 3], sa_z_1000, so_z_1000, A_1000, table_file, cloud[::2, ::2])
            atmos_ref_all_1000.append(atmos_ref2)

    atmos_ref_all_500 = np.array(atmos_ref_all_500).astype(float)
    atmos_ref_all_1000 = np.array(atmos_ref_all_1000).astype(float)

    blue = atmos_ref_all_500[0]
    red = atmos_ref_all_500[2]
    nir = atmos_ref_all_500[3]
    nir2 = cd.expand_band(atmos_ref_all_1000[0])

    f = gdal.Open('./ndata/IGBP_lai_500.tif')
    landcover = f.ReadAsArray()

    rvi = zs.RVI(nir, red)
    ndvi = zs.NDVI(nir, red)
    savi = zs.SAVI(nir, red, 0.5)
    gvi = zs.GVI(blue, red, nir, nir2)
    evi = zs.EVI(nir, red, blue)
    lai = zs.LAI(ndvi, landcover)

    for i in range(rvi.shape[0]):
        for j in range(rvi.shape[1]):
            if cloud[i][j] == 1:
                rvi[i][j] = -999
                ndvi[i][j] = -999
                savi[i][j] = -999
                gvi[i][j] = -999
                evi[i][j] = -999
                lai[i][j] = -999
    # 将地表生物量瞬时产品保存为nc文件
    out_head2 = '(FY3D-MERSI)'
    cd.zhishu_nc(rvi, lat_nc_500, lon_nc_500, out_dir, year, month, day, out_head2, name='RVI')
    cd.zhishu_nc(ndvi, lat_nc_500, lon_nc_500, out_dir, year, month, day, out_head2, name='NDVI')
    cd.zhishu_nc(savi, lat_nc_500, lon_nc_500, out_dir, year, month, day, out_head2, name='SAVI')
    cd.zhishu_nc(gvi, lat_nc_500, lon_nc_500, out_dir, year, month, day, out_head2, name='GVI')
    cd.zhishu_nc(evi, lat_nc_500, lon_nc_500, out_dir, year, month, day, out_head2, name='EVI')
    cd.zhishu_nc(lai, lat_nc_500, lon_nc_500, out_dir, year, month, day, out_head2, name='LAI')


def FY3D_MERSI_MonProducts(Ins_Dir, Mon_SNR, Mon_SAT, Mon_DT, Mon_Dir):
    file_list = os.listdir(Ins_Dir)
    file_list_rvi = list()
    file_list_ndvi = list()
    file_list_savi = list()
    file_list_gvi = list()
    file_list_evi = list()
    file_list_lai = list()

    for file in file_list:
        if file.find('RVI') != -1:
            file_list_rvi.append(file)
        if file.find('NDVI') != -1:
            file_list_ndvi.append(file)
        if file.find('SAVI') != -1:
            file_list_savi.append(file)
        if file.find('GVI') != -1:
            file_list_gvi.append(file)
        if file.find('EVI') != -1:
            file_list_evi.append(file)
        if file.find('LAI') != -1:
            file_list_lai.append(file)

    fn = Ins_Dir + '/' + file_list_rvi[0]

    data = nc.Dataset(fn, 'r')
    lat = data['RVI'].variables['lat'][:]
    lon = data['RVI'].variables['lon'][:]
    time = fn[-13:-3]
    time = time[:7]

    # 根据多天的地表生物量瞬时产品计算地表生物量月产品
    rvi = list()
    for f in file_list_rvi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['RVI'].variables['RVI'][:]
        rvi.append(data)
    rvi = np.array(rvi).astype(float)

    header = ['Mon-RVI(FY3D-MERSI)', 'RVI', time]
    rvi_aver = cd.zs_month(rvi, Mon_Dir, header, lat, lon)

    ndvi = list()
    for f in file_list_ndvi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['NDVI'].variables['NDVI'][:]
        ndvi.append(data)
    ndvi = np.array(ndvi).astype(float)

    header = ['Mon-NDVI(FY3D-MERSI)', 'NDVI', time]
    ndvi_aver = cd.zs_month(ndvi, Mon_Dir, header, lat, lon)

    savi = list()
    for f in file_list_savi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['SAVI'].variables['SAVI'][:]
        savi.append(data)
    savi = np.array(savi).astype(float)

    header = ['Mon-SAVI(FY3D-MERSI)', 'SAVI', time]
    savi_aver = cd.zs_month(savi, Mon_Dir, header, lat, lon)

    gvi = list()
    for f in file_list_gvi:
        fn = Ins_Dir + '/' + f

        data = nc.Dataset(fn, 'r')
        data = data['GVI'].variables['GVI'][:]
        gvi.append(data)
    gvi = np.array(savi).astype(float)

    header = ['Mon-GVI(FY3D-MERSI)', 'GVI', time]
    gvi_aver = cd.zs_month(gvi, Mon_Dir, header, lat, lon)

    evi = list()
    for f in file_list_evi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['EVI'].variables['EVI'][:]
        evi.append(data)
    evi = np.array(evi).astype(float)

    header = ['Mon-EVI(FY3D-MERSI)', 'EVI', time]
    evi_aver = cd.zs_month(evi, Mon_Dir, header, lat, lon)

    lai = list()
    for f in file_list_lai:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['LAI'].variables['LAI'][:]
        lai.append(data)
    lai = np.array(lai).astype(float)

    header = ['Mon-LAI(FY3D-MERSI)', 'LAI', time]
    lai_aver = cd.zs_month(lai, Mon_Dir, header, lat, lon)

    header = ['Mon-GPP(FY3D-MERSI)', 'Mon-NPP(FY3D-MERSI)']
    IGBP_file = './ndata/IGBP_2020.tif'

    mgnt.GPP_month(time[:4], time[5:], ndvi_aver, 7 / ndvi.shape[1], Mon_SAT, Mon_DT, Mon_SNR, IGBP_file, Mon_Dir,
                   header)


def FY3C_VIRR_MonProducts(Ins_Dir, Mon_SNR, Mon_SAT, Mon_DT, Mon_Dir):
    file_list = os.listdir(Ins_Dir)
    file_list_rvi = list()
    file_list_ndvi = list()
    file_list_savi = list()
    file_list_gvi = list()
    file_list_evi = list()
    file_list_lai = list()

    for file in file_list:
        if file.find('RVI') != -1:
            file_list_rvi.append(file)
        if file.find('NDVI') != -1:
            file_list_ndvi.append(file)
        if file.find('SAVI') != -1:
            file_list_savi.append(file)
        if file.find('GVI') != -1:
            file_list_gvi.append(file)
        if file.find('EVI') != -1:
            file_list_evi.append(file)
        if file.find('LAI') != -1:
            file_list_lai.append(file)

    fn = Ins_Dir + '/' + file_list_rvi[0]

    data = nc.Dataset(fn, 'r')
    lat = data['RVI'].variables['lat'][:]
    lon = data['RVI'].variables['lon'][:]
    time = fn[-13:-3]
    time = time[:7]

    # 根据多天的地表生物量瞬时产品计算地表生物量月产品
    rvi = list()
    for f in file_list_rvi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['RVI'].variables['RVI'][:]
        rvi.append(data)
    rvi = np.array(rvi).astype(float)

    header = ['Mon-RVI(FY3C-VIRR)', 'RVI', time]
    rvi_aver = cd.zs_month(rvi, Mon_Dir, header, lat, lon)

    ndvi = list()
    for f in file_list_ndvi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['NDVI'].variables['NDVI'][:]
        ndvi.append(data)
    ndvi = np.array(ndvi).astype(float)

    header = ['Mon-NDVI(FY3C-VIRR)', 'NDVI', time]
    ndvi_aver = cd.zs_month(ndvi, Mon_Dir, header, lat, lon)

    savi = list()
    for f in file_list_savi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['SAVI'].variables['SAVI'][:]
        savi.append(data)
    savi = np.array(savi).astype(float)

    header = ['Mon-SAVI(FY3C-VIRR)', 'SAVI', time]
    savi_aver = cd.zs_month(savi, Mon_Dir, header, lat, lon)

    gvi = list()
    for f in file_list_gvi:
        fn = Ins_Dir + '/' + f

        data = nc.Dataset(fn, 'r')
        data = data['GVI'].variables['GVI'][:]
        gvi.append(data)
    gvi = np.array(savi).astype(float)

    header = ['Mon-GVI(FY3C-VIRR)', 'GVI', time]
    gvi_aver = cd.zs_month(gvi, Mon_Dir, header, lat, lon)

    evi = list()
    for f in file_list_evi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['EVI'].variables['EVI'][:]
        evi.append(data)
    evi = np.array(evi).astype(float)

    header = ['Mon-EVI(FY3C-VIRR)', 'EVI', time]
    evi_aver = cd.zs_month(evi, Mon_Dir, header, lat, lon)

    lai = list()
    for f in file_list_lai:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['LAI'].variables['LAI'][:]
        lai.append(data)
    lai = np.array(lai).astype(float)

    header = ['Mon-LAI(FY3C-VIRR)', 'LAI', time]
    lai_aver = cd.zs_month(lai, Mon_Dir, header, lat, lon)

    header = ['Mon-GPP(FY3C-VIRR)', 'Mon-NPP(FY3C-VIRR)']
    IGBP_file = './ndata/IGBP_2020.tif'
    mgnt.GPP_month(time[:4], time[5:], ndvi_aver, 0.01, Mon_SAT, Mon_DT, Mon_SNR, IGBP_file, Mon_Dir, header)


def JPSS1_VIIRS_MonProducts(Ins_Dir, Mon_SNR, Mon_SAT, Mon_DT, Mon_Dir):
    file_list = os.listdir(Ins_Dir)
    file_list_rvi = list()
    file_list_ndvi = list()
    file_list_savi = list()
    file_list_gvi = list()
    file_list_evi = list()
    file_list_lai = list()

    for file in file_list:
        if file.find('RVI') != -1:
            file_list_rvi.append(file)
        if file.find('NDVI') != -1:
            file_list_ndvi.append(file)
        if file.find('SAVI') != -1:
            file_list_savi.append(file)
        if file.find('GVI') != -1:
            file_list_gvi.append(file)
        if file.find('EVI') != -1:
            file_list_evi.append(file)
        if file.find('LAI') != -1:
            file_list_lai.append(file)
    fn = Ins_Dir + '/' + file_list_rvi[0]

    data = nc.Dataset(fn, 'r')
    lat = data['RVI'].variables['lat'][:]
    lon = data['RVI'].variables['lon'][:]
    time = fn[-13:-3]
    time = time[:7]

    # 根据多天的地表生物量瞬时产品计算地表生物量月产品
    rvi = list()
    for f in file_list_rvi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['RVI'].variables['RVI'][:]
        rvi.append(data)
    rvi = np.array(rvi).astype(float)

    header = ['Mon-RVI(JPSS1-VIIRS)', 'RVI', time]
    rvi_aver = cd.zs_month(rvi, Mon_Dir, header, lat, lon)

    ndvi = list()
    for f in file_list_ndvi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['NDVI'].variables['NDVI'][:]
        ndvi.append(data)
    ndvi = np.array(ndvi).astype(float)

    header = ['Mon-NDVI(JPSS1-VIIRS)', 'NDVI', time]
    ndvi_aver = cd.zs_month(ndvi, Mon_Dir, header, lat, lon)

    savi = list()
    for f in file_list_savi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['SAVI'].variables['SAVI'][:]
        savi.append(data)
    savi = np.array(savi).astype(float)

    header = ['Mon-SAVI(JPSS1-VIIRS)', 'SAVI', time]
    savi_aver = cd.zs_month(savi, Mon_Dir, header, lat, lon)

    gvi = list()
    for f in file_list_gvi:
        fn = Ins_Dir + '/' + f

        data = nc.Dataset(fn, 'r')
        data = data['GVI'].variables['GVI'][:]
        gvi.append(data)
    gvi = np.array(savi).astype(float)

    header = ['Mon-GVI(JPSS1-VIIRS)', 'GVI', time]
    gvi_aver = cd.zs_month(gvi, Mon_Dir, header, lat, lon)

    evi = list()
    for f in file_list_evi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['EVI'].variables['EVI'][:]
        evi.append(data)
    evi = np.array(evi).astype(float)

    header = ['Mon-EVI(JPSS1-VIIRS)', 'EVI', time]
    evi_aver = cd.zs_month(evi, Mon_Dir, header, lat, lon)

    lai = list()
    for f in file_list_lai:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['LAI'].variables['LAI'][:]
        lai.append(data)
    lai = np.array(lai).astype(float)

    header = ['Mon-LAI(JPSS1-VIIRS)', 'LAI', time]
    lai_aver = cd.zs_month(lai, Mon_Dir, header, lat, lon)

    header = ['Mon-GPP(JPSS1-VIIRS)', 'Mon-NPP(JPSS1-VIIRS)']
    IGBP_file = './ndata/IGBP_2020.tif'
    mgnt.GPP_month(time[:4], time[5:], ndvi_aver, 0.01, Mon_SAT, Mon_DT, Mon_SNR, IGBP_file, Mon_Dir, header)


def NPP_VIIRS_MonProducts(Ins_Dir, Mon_SNR, Mon_SAT, Mon_DT, Mon_Dir):
    file_list = os.listdir(Ins_Dir)
    file_list_rvi = list()
    file_list_ndvi = list()
    file_list_savi = list()
    file_list_gvi = list()
    file_list_evi = list()
    file_list_lai = list()

    for file in file_list:
        if file.find('RVI') != -1:
            file_list_rvi.append(file)
        if file.find('NDVI') != -1:
            file_list_ndvi.append(file)
        if file.find('SAVI') != -1:
            file_list_savi.append(file)
        if file.find('GVI') != -1:
            file_list_gvi.append(file)
        if file.find('EVI') != -1:
            file_list_evi.append(file)
        if file.find('LAI') != -1:
            file_list_lai.append(file)

    fn = Ins_Dir + '/' + file_list_rvi[0]

    data = nc.Dataset(fn, 'r')
    lat = data['RVI'].variables['lat'][:]
    lon = data['RVI'].variables['lon'][:]
    lat = np.array(lat)
    lon = np.array(lon)
    time = fn[-13:-3]
    time = time[:7]

    # 根据多天的地表生物量瞬时产品计算地表生物量月产品
    rvi = list()
    for f in file_list_rvi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['RVI'].variables['RVI'][:]
        rvi.append(data)
    rvi = np.array(rvi).astype(float)

    header = ['Mon-RVI(NPP-VIIRS)', 'RVI', time]
    rvi_aver = cd.zs_month(rvi, Mon_Dir, header, lat, lon)

    ndvi = list()
    for f in file_list_ndvi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['NDVI'].variables['NDVI'][:]
        ndvi.append(data)
    ndvi = np.array(ndvi).astype(float)

    header = ['Mon-NDVI(NPP-VIIRS)', 'NDVI', time]
    ndvi_aver = cd.zs_month(ndvi, Mon_Dir, header, lat, lon)

    savi = list()
    for f in file_list_savi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['SAVI'].variables['SAVI'][:]
        savi.append(data)
    savi = np.array(savi).astype(float)

    header = ['Mon-SAVI(NPP-VIIRS)', 'SAVI', time]
    savi_aver = cd.zs_month(savi, Mon_Dir, header, lat, lon)

    gvi = list()
    for f in file_list_gvi:
        fn = Ins_Dir + '/' + f

        data = nc.Dataset(fn, 'r')
        data = data['GVI'].variables['GVI'][:]
        gvi.append(data)
    gvi = np.array(savi).astype(float)

    header = ['Mon-GVI(NPP-VIIRS)', 'GVI', time]
    gvi_aver = cd.zs_month(gvi, Mon_Dir, header, lat, lon)

    evi = list()
    for f in file_list_evi:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['EVI'].variables['EVI'][:]
        evi.append(data)
    evi = np.array(evi).astype(float)

    header = ['Mon-EVI(NPP-VIIRS)', 'EVI', time]
    evi_aver = cd.zs_month(evi, Mon_Dir, header, lat, lon)

    lai = list()
    for f in file_list_lai:
        fn = Ins_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['LAI'].variables['LAI'][:]
        lai.append(data)
    lai = np.array(lai).astype(float)

    header = ['Mon-LAI(NPP-VIIRS)', 'LAI', time]
    lai_aver = cd.zs_month(lai, Mon_Dir, header, lat, lon)

    header = ['Mon-GPP(NPP-VIIRS)', 'Mon-NPP(NPP-VIIRS)']
    IGBP_file = './ndata/IGBP_2020.tif'
    mgnt.GPP_month(time[:4], time[5:], ndvi_aver, 0.005, Mon_SAT, Mon_DT, Mon_SNR, IGBP_file, Mon_Dir, header)


def FY3D_MERSI_YearProducts(Mon_Dir, Year_Dir):
    file_list = os.listdir(Mon_Dir)
    file_list_rvi = list()
    file_list_ndvi = list()
    file_list_savi = list()
    file_list_gvi = list()
    file_list_evi = list()
    file_list_lai = list()
    file_list_gpp = list()
    file_list_npp = list()

    for file in file_list:
        if file.find('RVI') != -1:
            file_list_rvi.append(file)
        if file.find('NDVI') != -1:
            file_list_ndvi.append(file)
        if file.find('SAVI') != -1:
            file_list_savi.append(file)
        if file.find('GVI') != -1:
            file_list_gvi.append(file)
        if file.find('EVI') != -1:
            file_list_evi.append(file)
        if file.find('LAI') != -1:
            file_list_lai.append(file)
        if file.find('GPP') != -1:
            file_list_gpp.append(file)
        if file.find('NPP') != -1:
            file_list_npp.append(file)

    fn = Mon_Dir + '/' + file_list_rvi[0]

    data = nc.Dataset(fn, 'r')
    lat = data['RVI'].variables['lat'][:]
    lon = data['RVI'].variables['lon'][:]
    time = fn[-10:-3]
    time = time[:4]
    # 根据多月的地表生物量瞬时产品计算地表生物量年产品
    rvi = list()
    for f in file_list_rvi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['RVI'].variables['RVI'][:]
        rvi.append(data)
    rvi = np.array(rvi).astype(float)

    header = ['Mon-RVI(FY3D_MERSI)', 'RVI', time]
    rvi_aver = cd.zs_year(rvi, Year_Dir, header, lat, lon)

    ndvi = list()
    for f in file_list_ndvi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['NDVI'].variables['NDVI'][:]
        ndvi.append(data)
    ndvi = np.array(ndvi).astype(float)

    header = ['Mon-NDVI(FY3D_MERSI)', 'NDVI', time]
    ndvi_aver = cd.zs_year(ndvi, Year_Dir, header, lat, lon)

    savi = list()
    for f in file_list_savi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['SAVI'].variables['SAVI'][:]
        savi.append(data)
    savi = np.array(savi).astype(float)

    header = ['Mon-SAVI(FY3D_MERSI)', 'SAVI', time]
    savi_aver = cd.zs_year(savi, Year_Dir, header, lat, lon)

    gvi = list()
    for f in file_list_gvi:
        fn = Mon_Dir + '/' + f

        data = nc.Dataset(fn, 'r')
        data = data['GVI'].variables['GVI'][:]
        gvi.append(data)
    gvi = np.array(savi).astype(float)

    header = ['Mon-GVI(FY3D_MERSI)', 'GVI', time]
    gvi_aver = cd.zs_year(gvi, Year_Dir, header, lat, lon)

    evi = list()
    for f in file_list_evi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['EVI'].variables['EVI'][:]
        evi.append(data)
    evi = np.array(evi).astype(float)

    header = ['Mon-EVI(FY3D_MERSI)', 'EVI', time]
    evi_aver = cd.zs_year(evi, Year_Dir, header, lat, lon)

    lai = list()
    for f in file_list_lai:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['LAI'].variables['LAI'][:]
        lai.append(data)
    lai = np.array(lai).astype(float)

    header = ['Mon-LAI(FY3D_MERSI)', 'EVI', time]
    lai_aver = cd.zs_year(lai, Year_Dir, header, lat, lon)

    gpp = list()
    for f in file_list_gpp:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['GPP'].variables['GPP'][:]
        gpp.append(data)
    gpp = np.array(gpp).astype(float)

    header = ['Mon-GPP(FY3D_MERSI)', 'GPP', time]
    gpp_aver = cd.zs_year(gpp, Year_Dir, header, lat, lon)

    npp = list()
    for f in file_list_npp:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['NPP'].variables['NPP'][:]
        npp.append(data)
    npp = np.array(npp).astype(float)

    header = ['Mon-NPP(FY3D_MERSI)', 'NPP', time]
    npp_aver = cd.zs_year(npp, Year_Dir, header, lat, lon)


def FY3C_VIRR_YearProducts(Mon_Dir, Year_Dir):
    file_list = os.listdir(Mon_Dir)
    file_list_rvi = list()
    file_list_ndvi = list()
    file_list_savi = list()
    file_list_gvi = list()
    file_list_evi = list()
    file_list_lai = list()
    file_list_gpp = list()
    file_list_npp = list()

    for file in file_list:
        if file.find('RVI') != -1:
            file_list_rvi.append(file)
        if file.find('NDVI') != -1:
            file_list_ndvi.append(file)
        if file.find('SAVI') != -1:
            file_list_savi.append(file)
        if file.find('GVI') != -1:
            file_list_gvi.append(file)
        if file.find('EVI') != -1:
            file_list_evi.append(file)
        if file.find('LAI') != -1:
            file_list_lai.append(file)
        if file.find('GPP') != -1:
            file_list_gpp.append(file)
        if file.find('NPP') != -1:
            file_list_npp.append(file)

    fn = Mon_Dir + '/' + file_list_rvi[0]

    # 根据多月的地表生物量瞬时产品计算地表生物量年产品
    data = nc.Dataset(fn, 'r')
    lat = data['RVI'].variables['lat'][:]
    lon = data['RVI'].variables['lon'][:]
    time = fn[-10:-3]
    time = time[:4]

    rvi = list()
    for f in file_list_rvi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['RVI'].variables['RVI'][:]
        rvi.append(data)
    rvi = np.array(rvi).astype(float)

    header = ['Mon-RVI(FY3C_VIRR)', 'RVI', time]
    rvi_aver = cd.zs_year(rvi, Year_Dir, header, lat, lon)

    ndvi = list()
    for f in file_list_ndvi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['NDVI'].variables['NDVI'][:]
        ndvi.append(data)
    ndvi = np.array(ndvi).astype(float)

    header = ['Mon-NDVI(FY3C_VIRR)', 'NDVI', time]
    ndvi_aver = cd.zs_year(ndvi, Year_Dir, header, lat, lon)

    savi = list()
    for f in file_list_savi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['SAVI'].variables['SAVI'][:]
        savi.append(data)
    savi = np.array(savi).astype(float)

    header = ['Mon-SAVI(FY3C_VIRR)', 'SAVI', time]
    savi_aver = cd.zs_year(savi, Year_Dir, header, lat, lon)

    gvi = list()
    for f in file_list_gvi:
        fn = Mon_Dir + '/' + f

        data = nc.Dataset(fn, 'r')
        data = data['GVI'].variables['GVI'][:]
        gvi.append(data)
    gvi = np.array(savi).astype(float)

    header = ['Mon-GVI(FY3C_VIRR)', 'GVI', time]
    gvi_aver = cd.zs_year(gvi, Year_Dir, header, lat, lon)

    evi = list()
    for f in file_list_evi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['EVI'].variables['EVI'][:]
        evi.append(data)
    evi = np.array(evi).astype(float)

    header = ['Mon-EVI(FY3C_VIRR)', 'EVI', time]
    evi_aver = cd.zs_year(evi, Year_Dir, header, lat, lon)

    lai = list()
    for f in file_list_lai:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['LAI'].variables['LAI'][:]
        lai.append(data)
    lai = np.array(lai).astype(float)

    header = ['Mon-LAI(FY3C_VIRR)', 'EVI', time]
    lai_aver = cd.zs_year(lai, Year_Dir, header, lat, lon)

    gpp = list()
    for f in file_list_gpp:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['GPP'].variables['GPP'][:]
        gpp.append(data)
    gpp = np.array(gpp).astype(float)

    header = ['Mon-GPP(FY3C_VIRR)', 'GPP', time]
    gpp_aver = cd.zs_year(gpp, Year_Dir, header, lat, lon)

    npp = list()
    for f in file_list_npp:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['NPP'].variables['NPP'][:]
        npp.append(data)
    npp = np.array(npp).astype(float)

    header = ['Mon-NPP(FY3C_VIRR)', 'NPP', time]
    npp_aver = cd.zs_year(npp, Year_Dir, header, lat, lon)


def JPSS1_VIIRS_YearProducts(Mon_Dir, Year_Dir):
    file_list = os.listdir(Mon_Dir)
    file_list_rvi = list()
    file_list_ndvi = list()
    file_list_savi = list()
    file_list_gvi = list()
    file_list_evi = list()
    file_list_lai = list()
    file_list_gpp = list()
    file_list_npp = list()

    for file in file_list:
        if file.find('RVI') != -1:
            file_list_rvi.append(file)
        if file.find('NDVI') != -1:
            file_list_ndvi.append(file)
        if file.find('SAVI') != -1:
            file_list_savi.append(file)
        if file.find('GVI') != -1:
            file_list_gvi.append(file)
        if file.find('EVI') != -1:
            file_list_evi.append(file)
        if file.find('LAI') != -1:
            file_list_lai.append(file)
        if file.find('GPP') != -1:
            file_list_gpp.append(file)
        if file.find('NPP') != -1:
            file_list_npp.append(file)

    fn = Mon_Dir + '/' + file_list_rvi[0]

    data = nc.Dataset(fn, 'r')
    lat = data['RVI'].variables['lat'][:]
    lon = data['RVI'].variables['lon'][:]
    time = fn[-10:-3]
    time = time[:4]

    # 根据多月的地表生物量瞬时产品计算地表生物量年产品
    rvi = list()
    for f in file_list_rvi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['RVI'].variables['RVI'][:]
        rvi.append(data)
    rvi = np.array(rvi).astype(float)

    header = ['Mon-RVI(JPSS1_VIIRS)', 'RVI', time]
    rvi_aver = cd.zs_year(rvi, Year_Dir, header, lat, lon)

    ndvi = list()
    for f in file_list_ndvi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['NDVI'].variables['NDVI'][:]
        ndvi.append(data)
    ndvi = np.array(ndvi).astype(float)

    header = ['Mon-NDVI(JPSS1_VIIRS)', 'NDVI', time]
    ndvi_aver = cd.zs_year(ndvi, Year_Dir, header, lat, lon)

    savi = list()
    for f in file_list_savi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['SAVI'].variables['SAVI'][:]
        savi.append(data)
    savi = np.array(savi).astype(float)

    header = ['Mon-SAVI(JPSS1_VIIRS)', 'SAVI', time]
    savi_aver = cd.zs_year(savi, Year_Dir, header, lat, lon)

    gvi = list()
    for f in file_list_gvi:
        fn = Mon_Dir + '/' + f

        data = nc.Dataset(fn, 'r')
        data = data['GVI'].variables['GVI'][:]
        gvi.append(data)
    gvi = np.array(savi).astype(float)

    header = ['Mon-GVI(JPSS1_VIIRS)', 'GVI', time]
    gvi_aver = cd.zs_year(gvi, Year_Dir, header, lat, lon)

    evi = list()
    for f in file_list_evi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['EVI'].variables['EVI'][:]
        evi.append(data)
    evi = np.array(evi).astype(float)

    header = ['Mon-EVI(JPSS1_VIIRS)', 'EVI', time]
    evi_aver = cd.zs_year(evi, Year_Dir, header, lat, lon)

    lai = list()
    for f in file_list_lai:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['LAI'].variables['LAI'][:]
        lai.append(data)
    lai = np.array(lai).astype(float)

    header = ['Mon-LAI(JPSS1_VIIRS)', 'EVI', time]
    lai_aver = cd.zs_year(lai, Year_Dir, header, lat, lon)

    gpp = list()
    for f in file_list_gpp:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['GPP'].variables['GPP'][:]
        gpp.append(data)
    gpp = np.array(gpp).astype(float)

    header = ['Mon-GPP(JPSS1_VIIRS)', 'GPP', time]
    gpp_aver = cd.zs_year(gpp, Year_Dir, header, lat, lon)

    npp = list()
    for f in file_list_npp:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['NPP'].variables['NPP'][:]
        npp.append(data)
    npp = np.array(npp).astype(float)

    header = ['Mon-NPP(JPSS1_VIIRS)', 'NPP', time]
    npp_aver = cd.zs_year(npp, Year_Dir, header, lat, lon)


def NPP_VIIRS_YearProducts(Mon_Dir, Year_Dir):
    file_list = os.listdir(Mon_Dir)
    file_list_rvi = list()
    file_list_ndvi = list()
    file_list_savi = list()
    file_list_gvi = list()
    file_list_evi = list()
    file_list_lai = list()
    file_list_lst = list()
    file_list_gpp = list()
    file_list_npp = list()

    for file in file_list:
        if file.find('RVI') != -1:
            file_list_rvi.append(file)
        if file.find('NDVI') != -1:
            file_list_ndvi.append(file)
        if file.find('SAVI') != -1:
            file_list_savi.append(file)
        if file.find('GVI') != -1:
            file_list_gvi.append(file)
        if file.find('EVI') != -1:
            file_list_evi.append(file)
        if file.find('LAI') != -1:
            file_list_lai.append(file)
        if file.find('LST') != -1:
            file_list_lst.append(file)
        if file.find('GPP') != -1:
            file_list_gpp.append(file)
        if file.find('NPP(NPP-VIIRS)') != -1:
            file_list_npp.append(file)

    fn = Mon_Dir + '/' + file_list_rvi[0]

    data = nc.Dataset(fn, 'r')
    lat = data['RVI'].variables['lat'][:]
    lon = data['RVI'].variables['lon'][:]
    time = fn[-10:-3]
    time = time[:4]

    # 根据多月的地表生物量瞬时产品计算地表生物量年产品
    rvi = list()
    for f in file_list_rvi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['RVI'].variables['RVI'][:]
        rvi.append(data)
    rvi = np.array(rvi).astype(float)

    header = ['Mon-RVI(NPP_VIIRS)', 'RVI', time]
    rvi_aver = cd.zs_year(rvi, Year_Dir, header, lat, lon)

    ndvi = list()
    for f in file_list_ndvi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['NDVI'].variables['NDVI'][:]
        ndvi.append(data)
    ndvi = np.array(ndvi).astype(float)

    header = ['Mon-NDVI(NPP_VIIRS)', 'NDVI', time]
    ndvi_aver = cd.zs_year(ndvi, Year_Dir, header, lat, lon)

    savi = list()
    for f in file_list_savi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['SAVI'].variables['SAVI'][:]
        savi.append(data)
    savi = np.array(savi).astype(float)

    header = ['Mon-SAVI(NPP_VIIRS)', 'SAVI', time]
    savi_aver = cd.zs_year(savi, Year_Dir, header, lat, lon)

    gvi = list()
    for f in file_list_gvi:
        fn = Mon_Dir + '/' + f

        data = nc.Dataset(fn, 'r')
        data = data['GVI'].variables['GVI'][:]
        gvi.append(data)
    gvi = np.array(savi).astype(float)

    header = ['Mon-GVI(NPP_VIIRS)', 'GVI', time]
    gvi_aver = cd.zs_year(gvi, Year_Dir, header, lat, lon)

    evi = list()
    for f in file_list_evi:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['EVI'].variables['EVI'][:]
        evi.append(data)
    evi = np.array(evi).astype(float)

    header = ['Mon-EVI(NPP_VIIRS)', 'EVI', time]
    evi_aver = cd.zs_year(evi, Year_Dir, header, lat, lon)

    lai = list()
    for f in file_list_lai:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['LAI'].variables['LAI'][:]
        lai.append(data)
    lai = np.array(lai).astype(float)

    header = ['Mon-LAI(NPP_VIIRS)', 'EVI', time]
    lai_aver = cd.zs_year(lai, Year_Dir, header, lat, lon)

    gpp = list()
    for f in file_list_gpp:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        data = data['GPP'].variables['GPP'][:]
        gpp.append(data)
    gpp = np.array(gpp).astype(float)

    header = ['Mon-GPP(NPP_VIIRS)', 'GPP', time]
    gpp_aver = cd.zs_year(gpp, Year_Dir, header, lat, lon)

    npp = list()
    for f in file_list_npp:
        fn = Mon_Dir + '/' + f
        data = nc.Dataset(fn, 'r')
        print(fn)
        data = data['NPP'].variables['NPP'][:]
        npp.append(data)
    npp = np.array(npp).astype(float)

    header = ['Mon-NPP(NPP_VIIRS)', 'NPP', time]
    npp_aver = cd.zs_year(npp, Year_Dir, header, lat, lon)
