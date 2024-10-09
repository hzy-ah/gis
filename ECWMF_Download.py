#!/usr/bin/env python
# created by ZhangXing, OCT 16, 2022
import cdsapi
import os
import re
import calendar
from time import time
import datetime


def decreased_months(dt):
    if (dt.month > 2):
        month = dt.month - 1
        year = dt.year
    else:
        month = dt.month - 1 + 12
        year = dt.year - 1

    day = min(dt.day, calendar.monthrange(year, month)[1])

    return dt.replace(year=year, month=month, day=day)


def downloadonefile(yearStart, monStart, yearEnd, monEnd, dirs):
    years = range(yearStart, yearEnd + 1)
    yr = [str(i) for i in years]

    months = range(monStart, monEnd + 1)
    mo = [str(j).zfill(2) for j in months]

    for iyr in range(len(yr)):
        for imo in range(len(mo)):
            print('=======  date:  ' + yr[iyr] + '.' + mo[imo] + '  =======')

            outpathT2 = dirs + r"ECWMF/{0}{1:0>2d}/SAT".format(yr[iyr], int(mo[imo]))  # 存放气温的每月数据文件夹
            outpathDT = dirs + r"ECWMF/{0}{1:0>2d}/DT".format(yr[iyr], int(mo[imo]))  # 存放露点温度的每月数据文件夹
            outpathSSR = dirs + r"ECWMF/{0}{1:0>2d}/NSR".format(yr[iyr], int(mo[imo]))  # 存放净太阳辐射的每月数据文件夹
            # 判断文件夹是否存在
            if not os.path.exists(outpathT2):
                os.makedirs(outpathT2)

            # 保存的文件名
            outfilenameT2 = outpathT2 + r"/SAT_ERA5_{0}_{1:0>2d}.nc".format(yr[iyr], int(mo[imo]))
            outfilenameDT = outpathDT + r"/DT_ERA5_{0}_{1:0>2d}.nc".format(yr[iyr], int(mo[imo]))
            outfilenameSSR = outpathSSR + r"/NSR_ERA5_{0}_{1:0>2d}.nc".format(yr[iyr], int(mo[imo]))
            c = cdsapi.Client()

            if (os.path.isfile(outfilenameT2)):  # 如果存在文件名则返回
                print(outfilenameT2 + "已存在")
            else:
                print("开始下载：" + outfilenameT2)

                # 下载2m的气温

                c.retrieve(
                    'reanalysis-era5-single-levels-monthly-means',
                    {
                        'variable': '2m_temperature',  # set variables
                        'product_type': 'monthly_averaged_reanalysis',
                        'year': yr[iyr],
                        'month': mo[imo],
                        'time': [
                            '00:00'
                        ],
                        "area": [29, 115, 22, 122],  # set area, default: global
                        'format': 'netcdf'  # set output tpye, default: grib
                    },
                    outfilenameT2)  # set filename

            # 下载地面2M露点温度
            if not os.path.exists(outpathDT):
                os.makedirs(outpathDT)
            if (os.path.isfile(outfilenameDT)):  # 如果存在文件名则返回
                print(outfilenameDT + "已存在")
            else:
                print("开始下载：" + outfilenameDT)
                c.retrieve(
                    'reanalysis-era5-single-levels-monthly-means',
                    {
                        'variable': '2m_dewpoint_temperature',  # set variables
                        'product_type': 'monthly_averaged_reanalysis',
                        'year': yr[iyr],
                        'month': mo[imo],
                        'time': [
                            '00:00'
                        ],
                        "area": [29, 115, 22, 122],  # set area, default: global
                        'format': 'netcdf'  # set output tpye, default: grib
                    },
                    outfilenameDT)  # set filename
            # #下载地面太阳净辐射
            if not os.path.exists(outpathSSR):
                os.makedirs(outpathSSR)
            if (os.path.isfile(outpathSSR)):  # 如果存在文件名则返回
                print(outpathSSR + "已存在")
            else:
                print("开始下载：" + outfilenameSSR)
                c.retrieve(
                    'reanalysis-era5-single-levels-monthly-means',
                    {
                        'variable': 'surface_net_solar_radiation_clear_sky',  # set variables
                        'product_type': 'monthly_averaged_reanalysis',
                        'year': yr[iyr],
                        'month': mo[imo],
                        'time': [
                            '00:00'
                        ],
                        "area": [29, 115, 22, 122],  # set area, default: global
                        'format': 'netcdf'  # set output tpye, default: grib
                    },
                    outfilenameSSR)  # set filename
    print("运行完成!")


def ECWMF_ERA5_MetroDataDownload(Date_begin, Date_end, ECWMF_Dir):
    # 判断时间格式，不符合参数格式则程序退出！
    res_begin = re.search('^[0-9]{4}-[0-9]{2}-[0-9]{2}$', Date_begin)
    res_end = re.search('^[0-9]{4}-[0-9]{2}-[0-9]{2}$', Date_end)
    if res_begin is None or res_end is None:
        print('Please check your input date!')
        return False
    yearStart = int(Date_begin[0:4])
    monStart = int(Date_begin[5:7])
    yearEnd = int(Date_end[0:4])
    monEnd = int(Date_end[5:7])
    ts = time()
    begin = datetime.date(yearStart, monStart, 1)
    end = datetime.date(yearEnd, monEnd, 1)

    if begin <= end:
        downloadonefile(yearStart, monStart, yearEnd, monEnd, ECWMF_Dir)
    else:
        print("日期输入有误，请检查日期有效性！")
