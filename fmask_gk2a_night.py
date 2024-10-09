import numpy as np


def classification(data, SunZenith):  # 非夜间云做掩码处理
    SunZenith = SunZenith.astype(float)
    invalid_1 = np.logical_and(SunZenith < 75, SunZenith)
    SunZenith[invalid_1] = np.nan
    day = np.ma.masked_array(data, np.isnan(SunZenith))  # 小于75的太阳天顶角值做掩码处理（剩余为大于75的值）
    invalid_2 = np.logical_and(SunZenith > 75, SunZenith)
    SunZenith[invalid_2] = np.nan
    night = np.ma.masked_array(data, np.isnan(SunZenith))
    return day, night


def cloud_detection_gk2a_night(SunZenith, BT_7_cloudy_n, BT_14_cloudy_n,
                               BT_15_cloudy_n):  # 输入数据为  太阳高度角数据;7,14,15波段的亮温数据
    # 根据sunzenith太阳高度角进行筛选
    # BT_7_cloudy_d, BT_7_cloudy_n = classification(BT_7_cloudy, SunZenith)
    # BT_14_cloudy_d, BT_14_cloudy_n = classification(BT_14_cloudy, SunZenith)
    # BT_15_cloudy_d, BT_15_cloudy_n = classification(BT_15_cloudy, SunZenith)

    # 阈值调整
    threshold_1 = 0.5
    threshold_2 = 2.5
    threshold_3 = -6
    threshold_4 = 32

    condition_bt1 = BT_14_cloudy_n - BT_15_cloudy_n
    condition_bt2 = BT_7_cloudy_n - BT_14_cloudy_n > threshold_3
    condition_bt3 = BT_7_cloudy_n - BT_15_cloudy_n > threshold_4
    condition_bt_1 = np.logical_and(condition_bt1 > threshold_1, condition_bt1 < threshold_2)
    condition_bt_2 = np.logical_and(condition_bt2,
                                    (BT_7_cloudy_n - BT_14_cloudy_n) <= (BT_7_cloudy_n - BT_14_cloudy_n).max())
    condition_bt_3 = np.logical_and(condition_bt3,
                                    (BT_7_cloudy_n - BT_15_cloudy_n) <= (BT_7_cloudy_n - BT_15_cloudy_n).max())
    # 输出cloud,云为ture，非云为false
    cloud_gk2a_n = np.logical_or(np.logical_or(condition_bt_1, condition_bt_2), condition_bt_3)
    return cloud_gk2a_n


def Get_isCloud_Night(so_z, BT_7, BT_14, BT_15):
    isCloud = np.zeros((BT_7.shape[0], BT_7.shape[1]))
    for i in range(BT_7.shape[0]):
        for j in range(BT_7.shape[1]):
            isCloud[i][j] = cloud_detection_gk2a_night(so_z, BT_7[i][j], BT_14[i][j], BT_15[i][j])

    return isCloud
