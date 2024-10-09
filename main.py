import os
from Satellite_algorithm import *
from Satellite_algorithm_H import *
from ECWMF_Download import ECWMF_ERA5_MetroDataDownload

# L1B_Dir = '../data/AQUA'
# Date = '2023-04-02'
# Ins_Dir = './Data/Aqua'

# L1B_Dir = '../data/TERRA'
# Date = '2023-04-03'
# Ins_Dir = './Data/Terra'

L1B_Dir = '../data1/NPP/04'
Date = '2023-04-10'
Ins_Dir = './Data/JPSS1'

# L1B_Dir = '../data1'
# Date = '2023-04-16'
# Ins_Dir = '../data1/results'

# L1B_Dir = '../data/NOAA19'
# Date = '2023-04-10'
# Ins_Dir = './Data/NOAA19'

# L1B_Dir = '../Data/GK2A/00/05'
# DateTime = '2023-05-26 13:50'
# Ins_Dir = './Data/Ins_Gk-2A'
#
# L1B_Dir = '../data1/FY3C'
# Date = '2023-05-19'
# Ins_Dir = './Data/FY3C'

# L1B_Dir = '../data1/FY3D/04'
# Date = '2023-04-16'
# Ins_Dir = './Data/FY3D'

# L1B_Dir = 'F:/new/05/01'
# Date = '2023-04-16'
# Ins_Dir = 'F:/new/05/02'

# L1B_Dir = 'F:/new/06/01'
# Date = '2023-04-09'
# Ins_Dir = 'F:/new/06/02'

# Ins_Dir = './Data/Ins_Gk-2A/19'
# Day_Dir = './Data/Day_GK-2A'

# Date_begin = '2023-05-01'
# Date_end = '2023-05-30'
# ECWMF_Dir = './'

# Ins_Dir = './Data/Aqua'
# Mon_SAT = './ECWMF/202305/SAT/SAT_ERA5_2023_05.nc'
# Mon_DT = './ECWMF/202305/DT/DT_ERA5_2023_05.nc'
# Mon_NSR = './ECWMF/202305/NSR/NSR_ERA5_2023_05.nc'
# Mon_Dir = './Data/Mon_Aqua'

# Ins_Dir = '../Data3'
# Mon_SAT = './ECWMF/202305/SAT/SAT_ERA5_2023_05.nc'
# Mon_DT = './ECWMF/202305/DT/DT_ERA5_2023_05.nc'
# Mon_NSR = './ECWMF/202305/NSR/NSR_ERA5_2023_05.nc'
# Mon_Dir = '../Data3/Mon_Aqua'

# Ins_Dir = './Data/Terra'
# Mon_SAT = './ECWMF/202304/SAT/SAT_ERA5_2023_04.nc'
# Mon_DT = './ECWMF/202304/DT/DT_ERA5_2023_04.nc'
# Mon_NSR = './ECWMF/202304/NSR/NSR_ERA5_2023_04.nc'
# Mon_Dir = './Data/Mon_Terra'

# Ins_Dir = './Data/NOAA19/04'
# Mon_SAT = './ECWMF/202304/SAT/SAT_ERA5_2023_04.nc'
# Mon_DT = './ECWMF/202304/DT/DT_ERA5_2023_04.nc'
# Mon_NSR = './ECWMF/202304/NSR/NSR_ERA5_2023_04.nc'
# Mon_Dir = './Data/Mon_NOAA19'

# Ins_Dir = './Data/Day_GK-2A'
# Mon_SAT = './ECWMF/202305/SAT/SAT_ERA5_2023_05.nc'
# Mon_DT = './ECWMF/202305/DT/DT_ERA5_2023_05.nc'
# Mon_NSR = './ECWMF/202305/NSR/NSR_ERA5_2023_05.nc'
# Mon_Dir = './Data/Mon_GK-2A'

# Ins_Dir = './Data/FY3C'
# Mon_SAT = './ECWMF/202305/SAT/SAT_ERA5_2023_05.nc'
# Mon_DT = './ECWMF/202305/DT/DT_ERA5_2023_05.nc'
# Mon_NSR = './ECWMF/202305/NSR/NSR_ERA5_2023_05.nc'
# Mon_Dir = './Data/Mon_FY3C'

# Ins_Dir = './Data/FY3D/05'
# Mon_SAT = './ECWMF/202305/SAT/SAT_ERA5_2023_05.nc'
# Mon_DT = './ECWMF/202305/DT/DT_ERA5_2023_05.nc'
# Mon_NSR = './ECWMF/202305/NSR/NSR_ERA5_2023_05.nc'
# Mon_Dir = './Data/Mon_FY3D'

# Ins_Dir = './Data/JPSS1'
# Mon_SAT = './ECWMF/202305/SAT/SAT_ERA5_2023_05.nc'
# Mon_DT = './ECWMF/202305/DT/DT_ERA5_2023_05.nc'
# Mon_NSR = './ECWMF/202305/NSR/NSR_ERA5_2023_05.nc'
# Mon_Dir = './Data/Mon_JPSS1'
#
# Ins_Dir = 'F:/new/05/02'
# Mon_SAT = './ECWMF/202304/SAT/SAT_ERA5_2023_04.nc'
# Mon_DT = './ECWMF/202304/DT/DT_ERA5_2023_04.nc'
# Mon_NSR = './ECWMF/202304/NSR/NSR_ERA5_2023_04.nc'
# Mon_Dir = 'F:/new/05/03'

# Mon_Dir = './Data/Mon_Aqua'
# Year_Dir = './Data/Year_Aqua'

# Mon_Dir = './Data/Mon_Terra'
# Year_Dir='./Data/Year_Terra'

# Mon_Dir = './Data/Mon_NOAA19'
# Year_Dir='./Data/Year_NOAA19'

# Mon_Dir = './Data/Mon_GK-2A'
# Year_Dir='./Data/Year_GK-2A'

# Mon_Dir = './Data/Mon_FY3C'
# Year_Dir = './Data/Year_FY3C'

# Mon_Dir = './Data/Mon_FY3D'
# Year_Dir = './Data/Year_FY3D'
#
# Mon_Dir = './Data/Mon_JPSS1'
# Year_Dir = './Data/Year_JPSS1'

# Mon_Dir = './Data/Mon_NPP'
# Year_Dir = './Data/Year_NPP'


# 按间距中的绿色按钮以运行脚本。
if __name__ == '__main__':
    print('产品')
    # 生产各卫星的地表生物量瞬时产品
    # Aqua_MODIS_InsProducts(L1B_Dir, Date, Ins_Dir)
    # Terra_MODIS_InsProducts(L1B_Dir, Date, Ins_Dir)
    # NOAA19_AVHRR_InsProducts(L1B_Dir, Date, Ins_Dir)
    # GK2A_AMI_InsProducts(L1B_Dir, DateTime, Ins_Dir)
    # FY3C_VIRR_InsProducts(L1B_Dir, Date, Ins_Dir)
    # FY3D_MERSI_InsProducts(L1B_Dir, Date, Ins_Dir)
    JPSS1_VIIRS_InsProducts(L1B_Dir, Date, Ins_Dir)
    # NPP_VIIRS_InsProducts(L1B_Dir, Date, Ins_Dir)
    # GK2A_AMI_DayProducts(Ins_Dir, Day_Dir)

    # 通过ECWMF下载
    # ECWMF_ERA5_MetroDataDownload(Date_begin, Date_end, ECWMF_Dir)

    # 生产各卫星的地表生物量月产品
    # Aqua_MODIS_MonProducts(Ins_Dir, Mon_NSR, Mon_SAT, Mon_DT, Mon_Dir)
    # Terra_MODIS_MonProducts(Ins_Dir, Mon_NSR, Mon_SAT, Mon_DT, Mon_Dir)
    # NOAA19_AVHRR_MonProducts(Ins_Dir, Mon_NSR, Mon_SAT, Mon_DT, Mon_Dir)
    # GK2A_AMI_MonProducts(Ins_Dir, Mon_NSR, Mon_SAT, Mon_DT, Mon_Dir)
    # FY3C_VIRR_MonProducts(Ins_Dir, Mon_NSR, Mon_SAT, Mon_DT, Mon_Dir)
    # FY3D_MERSI_MonProducts(Ins_Dir, Mon_NSR, Mon_SAT, Mon_DT, Mon_Dir)
    # JPSS1_VIIRS_MonProducts(Ins_Dir, Mon_NSR, Mon_SAT, Mon_DT, Mon_Dir)
    # NPP_VIIRS_MonProducts(Ins_Dir, Mon_NSR, Mon_SAT, Mon_DT, Mon_Dir)

    # 生产各卫星的地表生物量年产品
    # Aqua_MODIS_YearProducts(Mon_Dir, Year_Dir)
    # Terra_MODIS_YearProducts(Mon_Dir, Year_Dir)
    # NOAA19_AVHRR_YearProducts(Mon_Dir, Year_Dir)
    # GK2A_AMI_YearProducts(Mon_Dir, Year_Dir)
    # FY3C_VIRR_YearProducts(Mon_Dir, Year_Dir)
    # FY3D_MERSI_YearProducts(Mon_Dir, Year_Dir)
    # JPSS1_VIIRS_YearProducts(Mon_Dir, Year_Dir)
    # NPP_VIIRS_YearProducts(Mon_Dir, Year_Dir)
