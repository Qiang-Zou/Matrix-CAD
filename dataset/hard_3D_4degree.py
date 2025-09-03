import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline

# Degree=5的knots和控制点
knots = [0,0,0,0,0, 
         0.0385,0.0769,0.1154,0.1538,0.1923,0.2308,0.2692,0.3077,0.3462,0.3846,
         0.4231,0.4615,0.5,0.5385,0.5769,0.6154,0.6538,0.6923,0.7308,0.7692,
         0.8077,0.8462,0.8846,0.9231,0.9615,
         1,1,1,1,1]  # 尾部节点 (p+1)

control_points = [
    [0.000000, 0.100000, 0.200000, 1.000000],
    [0.100000, 0.200000, 0.400000, 1.000000],
    [0.200000, 0.300000, 0.600000, 1.000000],
    [0.300000, 0.400000, 0.700000, 1.000000],
    [0.400000, 0.500000, 0.800000, 1.000000],
    [0.500000, 0.600000, 0.900000, 1.000000],
    [0.600000, 0.700000, 1.000000, 1.000000],
    [0.700000, 0.800000, 0.900000, 1.000000],
    [0.800000, 0.900000, 0.800000, 1.000000],
    [0.900000, 1.000000, 0.700000, 1.000000],
    [1.000000, 0.900000, 0.600000, 1.000000],
    [1.000000, 0.800000, 0.500000, 1.000000],
    [0.900000, 0.500000, 0.000000, 1.000000],
    [0.818437, 0.742070, 0.034483, 1.000000],
    [0.607011, 0.885420, 0.068966, 1.000000],
    [0.351945, 0.871591, 0.103448, 1.000000],
    [0.157257, 0.706222, 0.137931, 1.000000],
    [0.102345, 0.456752, 0.172414, 1.000000],
    [0.209602, 0.224920, 0.206897, 1.000000],
    [0.435287, 0.105269, 0.241379, 1.000000],
    [0.687363, 0.146595, 0.275862, 1.000000],
    [0.863030, 0.332044, 0.310345, 1.000000],
    [0.890648, 0.585988, 0.344828, 1.000000],
    [0.758955, 0.804865, 0.379310, 1.000000],
    [0.521656, 0.899413, 0.413793, 1.000000],
    [0.275525, 0.831076, 0.448276, 1.000000],
    [0.120939, 0.627721, 0.482759, 1.000000],
    [0.120939, 0.372279, 0.517241, 1.000000],
    [0.275525, 0.168924, 0.551724, 1.000000],
    [0.521656, 0.100587, 0.586207, 1.000000],
]

# Degree 5 B-Spline曲线
degree = 4
spl = BSpline(knots, control_points, degree)

# 在[0, 1]区间内生成100个样本点
x_vals = np.linspace(0, 1, 1000)
y_vals = spl(x_vals)

# 设置matplotlib风格为类似MATLAB的外观
# plt.style.use('fivethirtyeight')  # 使用 Matplotlib 自带的 ggplot 风格
plt.rcParams.update({
    'axes.linewidth': 1.5,           # 设置坐标轴线宽
    'axes.grid': True,               # 显示网格
    'grid.linestyle': '--',          # 网格线为虚线
    'grid.color': 'gray',            # 网格线颜色
    'axes.titlesize': 16,            # 标题字体大小
    'axes.labelsize': 14,            # 坐标轴标签字体大小
    'xtick.labelsize': 12,           # x轴刻度字体大小
    'ytick.labelsize': 12,           # y轴刻度字体大小
    'figure.figsize': (8, 6),        # 图像大小
    'font.family': 'serif',          # 字体样式
})

# 绘制 3D B-Spline 曲线
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 绘制 B-Spline 曲线
# ax.plot3D(y_vals[:, 0], y_vals[:, 1], y_vals[:, 2], label="B-Spline Degree=4", color='b', linewidth=2)
# ax.plot3D(y_vals[:, 0], y_vals[:, 1], y_vals[:, 2], label="B-Spline Degree=4", color=(6/255, 157/255, 255/255), linewidth=3)
# ax.plot3D(y_vals[:, 0], y_vals[:, 1], y_vals[:, 2], label="B-Spline Degree=4", color=(232/255, 210/255, 179/255), linewidth=3)
ax.plot3D(y_vals[:, 0], y_vals[:, 1], y_vals[:, 2], label="B-Spline Curve", color=(1/255, 84/255, 147/255), linewidth=3)
# 绘制控制点
# ax.scatter([p[0] for p in control_points], [p[1] for p in control_points], [p[2] for p in control_points], color='r', s=50, label="Control Points")
# ax.scatter([p[0] for p in control_points], [p[1] for p in control_points], [p[2] for p in control_points], color=(128/255, 128/255, 128/255), s=30, label="Control Points")
# ax.scatter([p[0] for p in control_points], [p[1] for p in control_points], [p[2] for p in control_points], color=(189/255, 154/255, 173/255), s=30, label="Control Points")
ax.scatter([p[0] for p in control_points], [p[1] for p in control_points], [p[2] for p in control_points], color=(201/255, 78/255, 101/255), s=30, label="Control Points",alpha=1)

# 绘制控制多边形
control_polygon_x = [p[0] for p in control_points]
control_polygon_y = [p[1] for p in control_points]
control_polygon_z = [p[2] for p in control_points]
# ax.plot3D(control_polygon_x, control_polygon_y, control_polygon_z, color='r', linestyle='--', linewidth=1.5, label="Control Polygon")
# ax.plot3D(control_polygon_x, control_polygon_y, control_polygon_z, color=(164/255, 224/255, 72/255), linestyle='-', linewidth=1.5, label="Control Polygon")
# ax.plot3D(control_polygon_x, control_polygon_y, control_polygon_z, color=(145/255, 147/255, 180/255), linestyle='-', linewidth=1.5, label="Control Polygon")
ax.plot3D(control_polygon_x, control_polygon_y, control_polygon_z, color=(201/255, 78/255, 101/255), linestyle='-', linewidth=1.5, label="Control Polygon")

# 设置标题和轴标签
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# # 设置标题字体大小
# ax.set_title("B-Spline Degree=4 in 3D", fontsize=16)

# 显示图例
# ax.legend(fontsize=20)

# 显示图形
plt.show()

# 保存为PNG图片
plt.savefig('hard_3D_4degree.pdf', bbox_inches='tight',pad_inches=0.2)
