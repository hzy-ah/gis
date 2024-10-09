from PIL import Image
import numpy as np
import netCDF4 as nc
from scipy import ndimage

from LST_h import *
import create_data as cd
import read_data as rd
import basic_function as bf
import h5py
import os
import zhishu as zs
import fmask_h
from osgeo import gdal
import math
def ec():
    # 假设原始29*29矩阵为image，这里使用随机数生成一个示例
    src_file1 = 'F:/Merge/ECWMF/202305/NSR/NSR_ERA5_2023_05.nc'
    src_file2 = 'F:/Merge/ECWMF/202305/DT/DT_ERA5_2023_05.nc'
    src_file3 = 'F:/Merge/ECWMF/202305/SAT/SAT_ERA5_2023_05.nc'
    Data_nc1 = nc.Dataset(src_file1)
    Data1 = Data_nc1['ssrc'][0][:]
    unit1 = Data_nc1['ssrc'].units
    image1 = np.array(Data1).astype(int)
    Data_nc2 = nc.Dataset(src_file2)
    Data2 = Data_nc2['d2m'][0][:]
    unit2 = Data_nc2['d2m'].units
    image2 = np.array(Data2).astype(int)
    Data_nc3 = nc.Dataset(src_file3)
    Data3 = Data_nc3['t2m'][0][:]
    unit3 = Data_nc3['t2m'].units
    image3 = np.array(Data3).astype(int)
    print(image1, unit1)

    # # 缩放后的图像尺寸
    # new_width, new_height = 700, 700
    #
    # image2 = ndimage.zoom(image, (700/29, 700/29), order=1)
    #
    # # 使用PIL（Pillow）创建新图像
    # pil_image = Image.fromarray(image)
    # new_pil_image = Image.fromarray(image2)
    #
    # # 保存原始图像为TIF格式
    # pil_image.save("F:/original_image.tif")
    #
    # # 保存扩充后的图像为TIF格式
    # new_pil_image.save("F:/resized_image.tif")

def xiaochu():
    print('1')

def remove_rows_with_value(matrix, value):
    # 创建一个布尔索引，表示每一行是否全是给定的value
    row_mask = np.all(matrix < value, axis=1)

    # 使用布尔索引删除满足条件的行
    new_matrix = matrix[~row_mask]

    # 返回删除的行的索引位置
    removed_rows = np.where(row_mask)[0]

    return new_matrix, removed_rows

def viirs():
    # geo_path = 'F:/new/05/GIMGO_npp_d20230416_t0455099_e0500503_b59415_c20230416052311498283_oeac_ops.h5'
    # data_path = 'F:/new/05/SVI01_npp_d20230416_t0455099_e0500503_b59415_c20230416052429818504_oeac_ops.h5'
    # data_path2 = 'F:/new/05/SVI02_npp_d20230416_t0455099_e0500503_b59415_c20230416052432089305_oeac_ops.h5'
    # data_path3 = 'F:/new/05/SVI03_npp_d20230416_t0455099_e0500503_b59415_c20230416052434443312_oeac_ops.h5'
    # out_dir = 'F:/new/05'

    # geo_path = 'F:/new/01/GIMGO_npp_d20230409_t0526375_e0532179_b59316_c20230409055402574797_oeac_ops.h5'
    # data_path = 'F:/new/01/SVI01_npp_d20230409_t0526375_e0532179_b59316_c20230409055517418239_oeac_ops.h5'
    # data_path2 = 'F:/new/01/SVI02_npp_d20230409_t0526375_e0532179_b59316_c20230409055519537722_oeac_ops.h5'
    # data_path3 = 'F:/new/01/SVI03_npp_d20230409_t0526375_e0532179_b59316_c20230409055521526044_oeac_ops.h5'
    # out_dir = 'F:/new/01'


    # geo_path = 'F:/new/07/GIMGO_j02_d20230409_t0458096_e0503502_b02127_c20230409052148464005_oeac_ops.h5'
    # data_path = 'F:/new/07/SVI01_j02_d20230409_t0458096_e0503502_b02127_c20230409052258652957_oeac_ops.h5'
    # data_path2 = 'F:/new/07/SVI02_j02_d20230409_t0458096_e0503502_b02127_c20230409052300720193_oeac_ops.h5'
    # data_path3 = 'F:/new/07/SVI03_j02_d20230409_t0458096_e0503502_b02127_c20230409052302697631_oeac_ops.h5'
    # out_dir = 'F:/new/07'

    geo_path = 'F:/new/03/GIMGO_j01_d20230416_t0546221_e0552021_b28019_c20230416060724644181_oeac_ops.h5'
    data_path = 'F:/new/03/SVI01_j01_d20230416_t0546221_e0552021_b28019_c20230416060841853989_oeac_ops.h5'
    data_path2 = 'F:/new/03/SVI02_j01_d20230416_t0546221_e0552021_b28019_c20230416060843784591_oeac_ops.h5'
    data_path3 = 'F:/new/03/SVI03_j01_d20230416_t0546221_e0552021_b28019_c20230416060845783933_oeac_ops.h5'
    out_dir = 'F:/new/03'

    # data_path = 'F:/data/JPSS/RNSCA-RVIRS_JPSS1_d20230501_t1223_svi01.h5'
    # data_path2 = 'F:/data/JPSS/RNSCA-RVIRS_JPSS1_d20230501_t1223_svi02.h5'
    # data_path3 = 'F:/data/JPSS/RNSCA-RVIRS_JPSS1_d20230501_t1223_svi03.h5'
    # geo_path = 'F:/data/JPSS/RNSCA-RVIRS_JPSS1_d20230501_t1223_gimgo.h5'
    # out_dir = 'F:/new/test2'
    h1, lat1, lon1, sa_z1, so_z1, A1 = bf.geo(geo_path, data_type='VIIRS')
    x_max = 0
    x_min = 10000
    y_max = 0
    y_min = 10000

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

    Ref_svi = []
    new_Ref_svi = []

    ref_svi01 = bf.band(data_path, data_type='VIIRS', band_type='ref_i', data_range=data_range)
    ref_svi02 = bf.band(data_path2, data_type='VIIRS', band_type='ref_i', data_range=data_range)
    ref_svi03 = bf.band(data_path3, data_type='VIIRS', band_type='ref_i', data_range=data_range)

    Ref_svi.append(ref_svi01)
    Ref_svi.append(ref_svi02)
    Ref_svi.append(ref_svi03)

    Ref_svi = np.array(Ref_svi).astype(float)


    # for i in range(Ref_svi.shape[1]):
    #     for j in range(Ref_svi.shape[2]):
    #         for num in range(Ref_svi.shape[0]):
    #             Ref_svi[num][i][j] = Ref_svi[num][i][j] / math.cos(so_z1[i][j] * math.pi / 180)


    value_to_remove = -800

    # 判断并删除包含给定值的行
    new_so_z1, removed_rows = remove_rows_with_value(so_z1, value_to_remove)
    for i in range(Ref_svi.shape[0]):
        new_ref = np.delete(Ref_svi[i], removed_rows, axis=0)
        new_Ref_svi.append(new_ref)


    new_Ref_svi = np.array(new_Ref_svi)

    new_lat1 = np.delete(lat1, removed_rows, axis=0)
    new_lon1 = np.delete(lon1, removed_rows, axis=0)
    new_sa_z1 = np.delete(sa_z1, removed_rows, axis=0)
    new_A1 = np.delete(A1, removed_rows, axis=0)

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

    for i in range(new_Ref_svi.shape[0]):
        ref_path1 = cd.create_tif(new_Ref_svi[i], out_dir, cd.varname(new_Ref_svi))
        cd.xml_file(ref_path1, lat_path1, lon_path1, vrt_path, new_Ref_svi[i].shape[0], new_Ref_svi[i].shape[1])
        bf.renaming(vrt_path)
        swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(new_Ref_svi), outputBounds, xRes, yRes)
        ref[i] = swath_data.ReadAsArray()
        ref[i] = cd.replace_fill_value(ref[i])
        swath_data = None
        os.remove(out)



    cd.create_tif(ref[0], out_dir, 'ref1')
    cd.create_tif(ref[1], out_dir, 'ref2')
    cd.create_tif(ref[2], out_dir, 'ref3')

    cd.xml_file(sa_z_path1, lat_path1, lon_path1, vrt_path, sa_z1.shape[0], sa_z1.shape[1])
    bf.renaming(vrt_path)
    swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(sa_z1), outputBounds, xRes, yRes)
    sa_z = swath_data.ReadAsArray()
    swath_data = None
    os.remove(out)

    cd.xml_file(so_z_path1, lat_path1, lon_path1, vrt_path, so_z1.shape[0], so_z1.shape[1])
    bf.renaming(vrt_path)
    swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(so_z1), outputBounds, xRes, yRes)
    so_z = swath_data.ReadAsArray()
    swath_data = None
    os.remove(out)

    cd.xml_file(A_path1, lat_path1, lon_path1, vrt_path, A1.shape[0], A1.shape[1])
    bf.renaming(vrt_path)
    swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(A1), outputBounds, xRes, yRes)
    A = swath_data.ReadAsArray()
    lat_nc, lon_nc = rd.read_lat_lon(swath_data)
    swath_data = None

    for i in range(sa_z.shape[0]):
        for j in range(sa_z.shape[1]):
            if sa_z[i][j] > 72:
                sa_z[i][j] = 0

            if so_z[i][j] > 80:
                so_z[i][j] = 0

            if A[i][j] > 180:
                A[i][j] = 0

    os.remove(out)
    red_cloud = ref[0]
    nir_cloud = ref[1]
    swir1_cloud = ref[2]
    ndvi_cloud = zs.NDVI(nir_cloud, red_cloud)
    cloud = fmask_h.fmask2(nir_cloud, red_cloud, swir1_cloud, ndvi_cloud)
    atmos_ref_all = list()
    for i in range(3):
        ref_i_whole = [1, 2, 3]
        ref_i_index = [0, 1, 2]
        table_file = './AT_LUTS/VIIRS_I' + str(ref_i_whole[i]) + '.nc'
        atmos_ref = bf.atmos(ref[int(ref_i_index[i])], sa_z, so_z, A, table_file, cloud)
        atmos_ref_all.append(atmos_ref)

    atmos_ref_all = np.array(atmos_ref_all).astype(float)
    cd.create_tif(atmos_ref_all[0], out_dir, 'srf1')
    cd.create_tif(atmos_ref_all[1], out_dir, 'srf2')
    cd.create_tif(atmos_ref_all[2], out_dir, 'srf3')

    red = atmos_ref_all[0]
    nir = atmos_ref_all[1]
    ndvi = zs.NDVI(nir, red)
    for i in range(ndvi.shape[0]):
        for j in range(ndvi.shape[1]):
            if cloud[i][j] == 1:
                ndvi[i][j] = 0
    cd.create_tif(ndvi, out_dir, 'ndvi')

    os.remove(lat_path1)
    os.remove(lon_path1)
    os.remove(sa_z_path1)
    os.remove(so_z_path1)
    os.remove(A_path1)


def tiff():
    path = 'F:/new/02/01'
    out_dir = 'F:/new/02/02'
    path_list = os.listdir(path)
    I1_all = []
    I2_all = []
    Lat = []
    Lon = []
    for L in path_list:
        file_path = path + '/' + L
        data = nc.Dataset(file_path)
        i1 = data.variables['375m Surface Reflectance Band I1'][:].data
        i2 = data.variables['375m Surface Reflectance Band I2'][:].data
        lat = data.variables['Latitude_at_375m_resolution'][:].data
        lon = data.variables['Longitude_at_375m_resolution'][:].data
        i1 = np.array(i1)
        i2 = np.array(i2)


        lat = np.array(lat)
        lon = np.array(lon)
        I1_all.append(i1)
        I2_all.append(i2)


        Lat.append(lat)
        Lon.append(lon)
    I1_all = np.array(I1_all)
    I2_all = np.array(I2_all)


    Lat = np.array(Lat)
    Lon = np.array(Lon)
    for i in range(I2_all.shape[0]-1):
        if i == 0:
            I1 = np.concatenate((I1_all[i][0:-20][:], I1_all[i + 1][:][:]), axis=0)
            I2 = np.concatenate((I2_all[i][0:-20][:], I2_all[i + 1][:][:]), axis=0)


            latitude = np.concatenate((Lat[i][0:-20][:], Lat[i + 1][:][:]), axis=0)
            longitude = np.concatenate((Lon[i][0:-20][:], Lon[i + 1][:][:]), axis=0)
        else:
            I1 = np.concatenate((I1[0:-20][:], I1_all[i + 1][:][:]), axis=0)
            I2 = np.concatenate((I2[0:-20][:], I2_all[i + 1][:][:]), axis=0)


            latitude = np.concatenate((latitude[0:-20][:], Lat[i + 1][:][:]), axis=0)
            longitude = np.concatenate((longitude[0:-20][:], Lon[i + 1][:][:]), axis=0)

    # print(I1.shape, ref.shape, longitude.shape, latitude.shape)

    x_max = 0
    x_min = 12000
    y_max = 0
    y_min = 12000

    for m in range(latitude.shape[0]):
        for n in range(latitude.shape[1]):
            if latitude[m][n] >= 20 and latitude[m][n] <= 30 and longitude[m][n] >= 114 and longitude[m][n] <= 124:
                if m > x_max:
                    x_max = m

                if m < x_min:
                    x_min = m

                if n > y_max:
                    y_max = n

                if n < y_min:
                    y_min = n

    latitude = latitude[x_min:x_max + 1, y_min:y_max + 1]
    longitude = longitude[x_min:x_max + 1, y_min:y_max + 1]
    I1 = I1[x_min:x_max + 1, y_min:y_max + 1]
    I2 = I2[x_min:x_max + 1, y_min:y_max + 1]

    ndvi = -9999 * np.ones((I1.shape[0], I1.shape[1]))
    for i in range(ndvi.shape[0]):
        for j in range(ndvi.shape[1]):
            if (I2[i][j] + I1[i][j]) == 0:
                ndvi[i][j] = -9999
            elif I2[i][j] != -9999 and I1[i][j] != -9999:
                ndvi[i][j] = (I2[i][j] - I1[i][j])/(I2[i][j] + I1[i][j])
    print(I2.shape)
    lat_path = cd.create_tif(latitude, out_dir, cd.varname(latitude))
    lon_path = cd.create_tif(longitude, out_dir, cd.varname(longitude))
    I2_path = cd.create_tif(I2, out_dir, cd.varname(I2))
    I1_path = cd.create_tif(I1, out_dir, cd.varname(I1))
    ndvi_path = cd.create_tif(ndvi, out_dir, cd.varname(ndvi))

    vrt_path = out_dir + '/' + 'vrt.xml'
    outputBounds = (115, 22, 122, 29)
    xRes = 0.005
    yRes = 0.005
    ref1 = np.zeros((1400, 1400))
    ref2 = np.zeros((1400, 1400))
    ndvi2 = np.zeros((1400, 1400))

    cd.xml_file(I1_path, lat_path, lon_path, vrt_path, I1.shape[0], I1.shape[1])
    bf.renaming(vrt_path)
    swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(I1), outputBounds, xRes, yRes)
    ref1 = swath_data.ReadAsArray()
    ref1 = cd.replace_fill_value(ref1)
    swath_data = None
    os.remove(out)
    cd.create_tif(ref1, out_dir, 'srf1')

    cd.xml_file(I2_path, lat_path, lon_path, vrt_path, I2.shape[0], I2.shape[1])
    bf.renaming(vrt_path)
    swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(I2), outputBounds, xRes, yRes)
    ref2 = swath_data.ReadAsArray()
    ref2 = cd.replace_fill_value(ref2)
    swath_data = None
    os.remove(out)
    cd.create_tif(ref2, out_dir, 'srf2')

    cd.xml_file(ndvi_path, lat_path, lon_path, vrt_path, ndvi.shape[0], ndvi.shape[1])
    bf.renaming(vrt_path)
    swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(ndvi), outputBounds, xRes, yRes)
    ndvi2 = swath_data.ReadAsArray()
    ndvi2 = cd.replace_fill_value(ndvi2)
    swath_data = None
    os.remove(out)
    cd.create_tif(ndvi2, out_dir, 'ndvi')

def LST_NPP():
    path = 'F:/new/04/03'
    out_dir = 'F:/new/04/04'
    path_list = os.listdir(path)
    Ref = []
    Lat = []
    Lon = []
    for L in path_list:
        file_path = path + '/' + L
        data = nc.Dataset(file_path)
        i2 = data.variables['VLST'][:].data
        lat = data.variables['Latitude'][:].data
        lon = data.variables['Longitude'][:].data
        i2 = np.array(i2)
        lat = np.array(lat)
        lon = np.array(lon)
        Ref.append(i2)
        Lat.append(lat)
        Lon.append(lon)
    Ref = np.array(Ref)
    Lat = np.array(Lat)
    Lon = np.array(Lon)
    for i in range(Ref.shape[0]-1):
        if i == 0:
            ref = np.concatenate((Ref[i][0:-2][:], Ref[i + 1][:][:]), axis=0)
            latitude = np.concatenate((Lat[i][0:-2][:], Lat[i + 1][:][:]), axis=0)
            longitude = np.concatenate((Lon[i][0:-2][:], Lon[i + 1][:][:]), axis=0)
        else:

            ref = np.concatenate((ref[0:-2][:], Ref[i + 1][:][:]), axis=0)
            latitude = np.concatenate((latitude[0:-2][:], Lat[i + 1][:][:]), axis=0)
            longitude = np.concatenate((longitude[0:-2][:], Lon[i + 1][:][:]), axis=0)


    print(ref.shape, longitude.shape, latitude.shape)

    x_max = 0
    x_min = 12000
    y_max = 0
    y_min = 12000

    for m in range(latitude.shape[0]):
        for n in range(latitude.shape[1]):
            if latitude[m][n] >= 20 and latitude[m][n] <= 30 and longitude[m][n] >= 114 and longitude[m][n] <= 124:
                if m > x_max:
                    x_max = m

                if m < x_min:
                    x_min = m

                if n > y_max:
                    y_max = n

                if n < y_min:
                    y_min = n

    latitude = latitude[x_min:x_max + 1, y_min:y_max + 1]
    longitude = longitude[x_min:x_max + 1, y_min:y_max + 1]
    ref = ref[x_min:x_max + 1, y_min:y_max + 1]

    for i in range(ref.shape[0]):
        for j in range(ref.shape[1]):
            if ref[i][j] != -32768:
                ref[i][j] = ref[i][j]*0.005 + 200



    lat_path = cd.create_tif(latitude, out_dir, cd.varname(latitude))
    lon_path = cd.create_tif(longitude, out_dir, cd.varname(longitude))
    ref_path = cd.create_tif(ref, out_dir, cd.varname(ref))

    vrt_path = out_dir + '/' + 'vrt.xml'
    outputBounds = (115, 22, 122, 29)
    xRes = 0.005
    yRes = 0.005
    ref2 = np.zeros((1400, 1400))




    cd.xml_file(ref_path, lat_path, lon_path, vrt_path, ref.shape[0], ref.shape[1])
    bf.renaming(vrt_path)
    swath_data, out = bf.swath_georeference_vrt(out_dir, cd.varname(ref), outputBounds, xRes, yRes)
    ref2 = swath_data.ReadAsArray()
    ref2 = cd.replace_fill_value(ref2)
    swath_data = None
    os.remove(out)
    cd.create_tif(ref2, out_dir, 'lst')

def quzao():
    out_dir = 'F:/new/07'
    dataset = gdal.Open('F:/new/07/lst-2.tif')
    band = dataset.GetRasterBand(1)
    data = band.ReadAsArray()
    IGBP_data = gdal.Open('F:/Merge/ndata/IGBP_LST.tif').ReadAsArray()
    # data = cd.replace_fill_value(data)
    # data = cd.replace_fill_value(data)
    # data = cd.replace_fill_value(data)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            data[i][j] = data[i][j] - 1.5
            if IGBP_data[i][j] == 17 or data[i][j] < 280:
                data[i][j] = -999
    cd.create_tif(data, out_dir, 'lst-21')



if __name__ == "__main__":
    # print('hello world!')
    # viirs()
    # tiff()
    # LST_NPP()
    quzao()
