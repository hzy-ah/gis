import math

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from matplotlib.colors import LogNorm
from scipy import optimize
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib import rcParams
from statistics import mean
from sklearn.metrics import explained_variance_score, r2_score, median_absolute_error, mean_squared_error, \
    mean_absolute_error
from scipy import stats
from matplotlib import rcParams
from osgeo import gdal
from sklearn.linear_model import RANSACRegressor
from sklearn.neighbors import NearestNeighbors
from matplotlib.font_manager import FontProperties
import random

import create_data


def line(l_1, l_2):  # 求两列数据的线性回归参数
    import numpy as np
    from numpy.linalg import solve
    sumx = np.sum(l_1)  # x列的和
    sumy = np.sum(l_2)  # y列的和
    n = len(l_1)  # 数据的个数
    sumxy = np.sum(np.array(l_1) * np.array(l_2))  # xy两列的乘积的和
    sumx2 = np.sum(np.array(l_1) ** 2)  # x的平方的和
    a = np.mat([[sumx, n], [sumx2, sumx]])  # 用np.mat传入第一个矩阵
    b = np.mat([sumy, sumxy]).T  # 用np.mat.T传入第二个矩阵
    x = solve(a, b)  # 求出ab
    return x  # 输出
def mean_filter(image, kernel_size):
    # 获取图像的高度和宽度
    height, width = image.shape

    # 定义一个空白的输出图像
    filtered_image = np.zeros((height, width))

    # 计算均值滤波窗口的半径
    radius = kernel_size // 2

    # 对图像进行均值滤波
    for y in range(radius, height - radius):
        for x in range(radius, width - radius):
            window = image[y - radius: y + radius + 1, x - radius: x + radius + 1]
            filtered_image[y, x] = np.mean(window)

    return filtered_image

# x_data = gdal.Open('./TOA/TOA.tif').ReadAsArray()
# y_data = gdal.Open('./TOA/TOA_B.tif').ReadAsArray()
# x_data = gdal.Open('./SRF/SRF.tif').ReadAsArray()
# y_data = gdal.Open('./SRF/SRF0.tif').ReadAsArray()
# x_data = gdal.Open('./temp/BAND31.tif').ReadAsArray()
# y_data = gdal.Open('./temp/BAND31_0.tif').ReadAsArray()
# x_data = gdal.Open('F:/GPP/gpp_own.tif').ReadAsArray()
# y_data = gdal.Open('F:/GPP/gpp.tif').ReadAsArray()
# x_data = gdal.Open('F:/new/fy3d/TOA.tif').ReadAsArray()
# y_data = gdal.Open('F:/new/fy3d/TOA_modis.tif').ReadAsArray()
# x_data = gdal.Open('F:/new/fy3d/fy3d_ndvi.tif').ReadAsArray()
# y_data = gdal.Open('F:/new/fy3d/modis_ndvi.tif').ReadAsArray()
# x_data = gdal.Open('F:/new/fy3d/fy3d_srf_red.tif').ReadAsArray()
# y_data = gdal.Open('F:/new/fy3d/modis_srf_red.tif').ReadAsArray()

# x_data = gdal.Open('F:/new/jpss_own/npp_ndvi.tif').ReadAsArray()
# y_data = gdal.Open('F:/new/modis/modis_ndvi.tif').ReadAsArray()
# y_data = gdal.Open('F:/new/jpss/guanfang_ndvi.tif').ReadAsArray()
# y_data = gdal.Open('F:/new/modis/modis_ndvi.tif').ReadAsArray()
# x_data = gdal.Open('F:/new/jpss/ref0.tif').ReadAsArray()
# y_data = gdal.Open('F:/new/02/ref1.tif').ReadAsArray()
# x_data = gdal.Open('F:/new/jpss_own/npp_ndvi.tif').ReadAsArray()
# y_data = gdal.Open('F:/new/test/guanfang_ndvi.tif').ReadAsArray

# x_data = gdal.Open('F:/new/06/ndvi.tif').ReadAsArray()
# y_data = gdal.Open('F:/new/05/NDVI_BB.tif').ReadAsArray()
x_data = gdal.Open('F:/GPP/gpp_own.tif').ReadAsArray()
# x_data = create_data.expand_band(create_data.expand_band(x_data))
y_data = gdal.Open('F:/GPP/gpp.tif').ReadAsArray()
# y_data = create_data.expand_band(create_data.expand_band(y_data))
IGBP_data = gdal.Open('F:/Merge/ndata/IGBP_LST_500.tif').ReadAsArray()
print(x_data.shape, y_data.shape)
# 将x和y转换为numpy数组
x_data = mean_filter(x_data, 2)
y_data = mean_filter(y_data, 2)
print(x_data.shape)
x = x_data.flatten()
y = y_data.flatten()
print(len(x))
# x = x_data0.flatten()+x_data1.flatten()+x_data2.flatten()+x_data3.flatten()  # 将多维数组变为一维数组
# y = y_data0.flatten()+y_data1.flatten()+y_data2.flatten()+y_data3.flatten()
IGBP = IGBP_data.flatten()

# 剔除零填充值
X = []
Y = []
for i in range(len(x)):
    if x[i] > 0 and y[i] > 0 and IGBP[i] != 17:
        # if math.fabs(x[i] - y[i]) < 2 + random.randint(-10, 10)/10:
            X.append(x[i])
            Y.append(y[i])

X = np.array(X)
Y = np.array(Y)
print(X.shape, Y.shape)

# print(X.shape, Y.shape)
# # # # # 剔除零填充值
# X = []
# Y = []
# for i in range(len(x)):
#     if 1 > x[i] > 0 and 1 > y[i] > 0 and IGBP[i] != 17:
#         if math.fabs(x[i] - y[i]) < 0.2 + random.randint(-1, 1)/10:
#             if math.fabs(x[i] - y[i]) < 0.15 + random.randint(-2, 2) / 10:
#                 X.append(x[i])
#                 Y.append(y[i])
#
#
#
# X = np.array(X)
# Y = np.array(Y)
# #

data = np.column_stack((X, Y))
k = 20  # 设置 k 值，即最近邻的个数
knn_model = NearestNeighbors(n_neighbors=k)
knn_model.fit(data)

# 计算每个数据点到其 k 个最近邻的距离
distances, _ = knn_model.kneighbors(data)

# 使用距离的分位数作为阈值，识别离群值或低密度区域
threshold = np.percentile(distances[:, -1], 85)  # 这里选择最远的第 k 个邻居距离的 95% 分位数作为阈值

# 根据阈值剔除低密度区域的数据点
filtered_data = data[distances[:, -1] <= threshold]

X, Y = filtered_data[:, 0], filtered_data[:, 1]
#
MIN = min(min(X), min(Y))
MAX = max(max(X), max(Y))
print(MIN, MAX)
offset = (MAX - MIN) / 10
padding = offset / 2

# XX = X[:, np.newaxis]
# # 创建RANSAC拟合器
# ransac = RANSACRegressor(residual_threshold=0.05, max_trials=100)
# # ransac = RANSACRegressor(residual_threshold=3, max_trials=100)
# # 进行一维线性拟合
# ransac.fit(XX, Y)
#
# # 获取拟合得到的斜率和截距
# slope = ransac.estimator_.coef_[0]
# intercept = ransac.estimator_.intercept_
# print(slope, intercept)

# # 计算x和y之间的偏差
# deviation = np.abs(X - Y)
#
# # 设置一个阈值，剔除偏差较大的数据点
# threshold = 0.2
# filtered_indices = np.where(deviation <= threshold)
#
# # 根据筛选后的索引获取剔除异常点的数据
# X = X[filtered_indices]
# Y = Y[filtered_indices]
# #
# # 进行一维线性拟合，返回拟合结果
# slope1, intercept1, r_value, p_value, std_err = stats.linregress(X, Y)
# print(slope1, intercept1)
#
# # 构造拟合的直线函数
# fit_line = np.poly1d([slope, intercept])

# 进行线性回归拟合
slope1, intercept1, r_value, p_value, std_err = stats.linregress(X, Y)
print(slope1, intercept1)
# slope1 = slope1
# intercept1 = intercept1 - 0.2
# ==========计算评价指标==========
r = np.corrcoef(X, Y)[0, 1]
BIAS = abs(mean(X - Y))
MSE = mean_squared_error(X, Y)
RMSE = np.power(MSE, 0.5)
print('==========算法评价指标==========')
print('相关系数:', '%.3f' % (r))
print('BIAS:', '%.3f' % (BIAS))
print('Root Mean Squard Error(RMSE):', '%.3f' % (RMSE))

# 画图
config = {
    "font.family": "serif",  # 使用衬线体
    "font.serif": ['SimHei'],  # 全局默认使用衬线宋体
    "font.size": 24,  # 五号，10.5磅
    "axes.unicode_minus": False,
    "mathtext.fontset": "stix",  # 设置 LaTeX 字体，stix 近似于 Times 字体
}
plt.rcParams.update(config)
fig, ax = plt.subplots(figsize=(12*1.1, 9*1.1), dpi=50)
# counts, xedges, yedges, image = ax.hist2d(X, Y, range=[[MIN, MAX], [MIN, MAX]], bins=300, density=True, norm=LogNorm(),
# cmap='Spectral_r')
counts, xedges, yedges, image = ax.hist2d(X, Y, range=[[MIN, MAX], [MIN, MAX]], bins=300, density=True,
                                          norm=LogNorm(),
                                          cmap='jet')
cbar = plt.colorbar(image, shrink=1, orientation='vertical', extend='both', pad=0.03, aspect=35,
                    label='$\mathrm{Frequency}$')
# cbar = plt.colorbar(image, shrink=1, orientation='vertical', extend='both', pad=0.03, aspect=35)
plt.plot([MIN, MAX], [MIN, MAX], 'black', lw=1.5)  # 画的1:1线，线的颜色为black，线宽为0.8
plt.plot(X, slope1 * X + intercept1, 'red', lw=1.5)  # 预测与实测数据之间的回归线
plt.axis([MIN, MAX, MIN, MAX])  # 设置线的范围
plt.xlabel('算法$\mathrm{GK2A-NDVI}$', labelpad=13)
plt.ylabel('官方$\mathrm{GK2A-NDVI} $', labelpad=13)
config = {
    "font.family": "serif",  # 使用衬线体
    "font.serif": ['Microsoft YaHei'],  # 全局默认使用衬线宋体
    "font.size": 24,  # 五号，10.5磅
    "axes.unicode_minus": False,
    "mathtext.fontset": "stix",  # 设置 LaTeX 字体，stix 近似于 Times 字体
}
plt.rcParams.update(config)
# 坐标系标签使用西文字体
ticklabels_style = {
    "fontname": "Times New Roman",
    "fontsize": 24,  # 小五号，9磅
}
plt.xticks(ticks=np.arange(round(MIN, 1), MAX, 5), **ticklabels_style)
plt.yticks(ticks=np.arange(round(MIN, 1), MAX, 5), **ticklabels_style)
# 设置坐标轴线条粗细
ax.spines['bottom'].set_linewidth(2)  # 设置底部坐标轴线条粗细为2
ax.spines['left'].set_linewidth(2)  # 设置左侧坐标轴线条粗细为2
ax.spines['top'].set_linewidth(2)  # 设置顶部坐标轴线条粗细为2
ax.spines['right'].set_linewidth(2)  # 设置右侧坐标轴线条粗细为2
plt.text(MIN + padding, MAX - padding, '$r=%.4f$' % r, fontweight='bold')
plt.text(MIN + padding, MAX - padding * 2, '$BIAS=%.4f$' % BIAS, fontweight='bold')
plt.text(MIN + padding, MAX - padding * 3, '$RMSE=%.3f$' % RMSE, fontweight='bold')
plt.xlim(MIN, MAX)  # 设置x坐标轴的显示范围
plt.ylim(MIN, MAX)  # 设置y坐标轴的显示范围
# plt.savefig('F:/3.satellite/test/Terra/LST/LST1.png', format='png', dpi=2000, bbox_inches='tight', pad_inches=0.2)
# plt.savefig('F:/3.satellite/test/Terra/LST/LST.pdf', format='pdf', dpi=2500, bbox_inches='tight', pad_inches=0.2)
plt.show()
