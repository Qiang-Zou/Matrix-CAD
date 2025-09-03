import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline

# Degree=5的knots和控制点
knots = [0,0,0,0,0,0,
        0.0286,0.0571,0.0857,0.1143,0.1429,0.1714,0.2,0.2286,0.2571,0.2857,
        0.3143,0.3429,0.3714,0.4,0.4286,0.4571,0.4857,0.5143,0.5429,0.5714,
        0.6,0.6286,0.6571,0.6857,0.7143,0.7429,0.7714,0.8,0.8286,0.8571,
        0.8857,0.9143,0.9429,0.9714,
        1,1,1,1,1,1]  # 尾部节点 (p+1)

control_points = [
[0.800000, 0.500000, 0.500000, 1.000000],
[0.707817, 0.716361, 0.540233, 1.000000],
[0.487920, 0.799757, 0.580206, 1.000000],
[0.275447, 0.698937, 0.619658, 1.000000],
[0.200973, 0.475860, 0.658334, 1.000000],
[0.310266, 0.267619, 0.695983, 1.000000],
[0.536161, 0.202187, 0.732362, 1.000000],
[0.739833, 0.319777, 0.767233, 1.000000],
[0.796115, 0.548123, 0.800371, 1.000000],
[0.670419, 0.746895, 0.831561, 1.000000],
[0.439992, 0.793937, 0.860601, 1.000000],
[0.246443, 0.660340, 0.887302, 1.000000],
[0.208717, 0.428205, 0.911492, 1.000000],
[0.350000, 0.240192, 0.933013, 1.000000],
[0.583465, 0.211845, 0.951725, 1.000000],
[0.765637, 0.360583, 0.967508, 1.000000],
[0.784561, 0.595000, 0.980259, 1.000000],
[0.628608, 0.771035, 0.989895, 1.000000],
[0.393619, 0.780505, 0.996354, 1.000000],
[0.224006, 0.617590, 0.999594, 1.000000],
[0.224006, 0.382410, 0.999594, 1.000000],
[0.393619, 0.219495, 0.996354, 1.000000],
[0.628608, 0.228965, 0.989895, 1.000000],
[0.784561, 0.405000, 0.980259, 1.000000],
[0.765637, 0.639417, 0.967508, 1.000000],
[0.583465, 0.788155, 0.951725, 1.000000],
[0.350000, 0.759808, 0.933013, 1.000000],
[0.208717, 0.571795, 0.911492, 1.000000],
[0.246443, 0.339660, 0.887302, 1.000000],
[0.439992, 0.206063, 0.860601, 1.000000],
[0.670419, 0.253105, 0.831561, 1.000000],
[0.796115, 0.451877, 0.800371, 1.000000],
[0.739833, 0.680223, 0.767233, 1.000000],
[0.536161, 0.797813, 0.732362, 1.000000],
[0.310266, 0.732381, 0.695983, 1.000000],
[0.200973, 0.524140, 0.658334, 1.000000],
[0.275447, 0.301063, 0.619658, 1.000000],
[0.487920, 0.200243, 0.580206, 1.000000],
[0.707817, 0.283639, 0.540233, 1.000000],
[0.800000, 0.500000, 0.500000, 1.000000]
]

print("Length of knots:", len(knots))
print("Shape of control_points:", np.array(control_points).shape)

# Degree 5 B-Spline曲线
degree = 5
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

plt.savefig('hard_3D_5degree.pdf', bbox_inches='tight',pad_inches=0.2)