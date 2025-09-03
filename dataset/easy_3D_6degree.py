import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline

# Degree=5的knots和控制点
knots = [0, 0, 0, 0, 0, 0, 0,   # 前导节点 (p+1)
        #  0.15, 0.3, 
        # #  0.45, 0.6, 0.75, 
        #  0.9,   # 内部节点
         1, 1, 1, 1, 1, 1,1]  # 尾部节点 (p+1)

control_points = [
    # [0.0, 0.1, 0.2, 1.0],
    # [0.1, 0.3, 0.4, 1.0],
    [0.2, 0.5, 0.6, 1.0],
    # [0.3, 0.7, 0.5, 1.0],
    [0.4, 0.9, 0.4, 1.0],
    [0.5, 1.0, 0.3, 1.0],
    # [0.6, 0.8, 0.2, 1.0],
    [0.7, 0.8, 0.3, 1.0],
    [0.8, 0.5, 0.4, 1.0],
    [0.9, 0.2, 0.5, 1.0],
    [1.0, 0.9, 0.6, 1.0],
    # [1.0, 0.1, 0.7, 1.0]
]

# Degree 5 B-Spline曲线
degree = 6
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
ax.plot3D(y_vals[:, 0], y_vals[:, 1], y_vals[:, 2], label="B-Spline Curve", color=(1/255, 84/255, 147/255), linewidth=3)
# 绘制控制点
ax.scatter([p[0] for p in control_points], [p[1] for p in control_points], [p[2] for p in control_points], color=(201/255, 78/255, 101/255), s=30, label="Control Points",alpha=1)

# 绘制控制多边形
control_polygon_x = [p[0] for p in control_points]
control_polygon_y = [p[1] for p in control_points]
control_polygon_z = [p[2] for p in control_points]
ax.plot3D(control_polygon_x, control_polygon_y, control_polygon_z, color=(201/255, 78/255, 101/255), linestyle='-', linewidth=1.5, label="Control Polygon")

# ax.set_title("B-Spline Degree=5 in 3D")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
# ax.legend(fontsize=20)
plt.show()

plt.savefig('easy_3D_6degree.pdf',bbox_inches='tight',pad_inches=0.2)