import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline

# Degree=5的knots和控制点
knots = [0, 0, 0, 0,0,   # 前导节点 (p+1)
         0.1, 0.2, 0.4, 0.5,   # 内部节点
         1, 1, 1, 1,1]  # 尾部节点 (p+1)

control_points = [
    [0.0, 0.0, 1.0, 1.0],
    [0.0, 0.5, 1.0, 1.0],
    [0.25, 0.5, 1.0, 1.0],
    [0.25, 1.0, 1.0, 1.0],
    [0.5, 1.0, 1.0, 1.0],
    [0.5, 0.5, 1.0, 1.0],
    [0.75, 0.5, 1.0, 1.0],
    [0.75, 0.0, 1.0, 1.0],
    [1.0, 0.0, 1.0, 1.0]
]

# Degree 5 B-Spline曲线
degree = 4
spl = BSpline(knots, control_points, degree)

# 在[0, 1]区间内生成100个样本点
x_vals = np.linspace(0, 1, 1000)
y_vals = spl(x_vals)

# 绘制2D B-Spline曲线
plt.figure(figsize=(8, 6))

# 曲线
plt.plot(y_vals[:, 0], y_vals[:, 1], label="B-Spline Curve",color=(1/255, 84/255, 147/255), linewidth=3)
# 控制点
plt.scatter([p[0] for p in control_points], [p[1] for p in control_points], color=(201/255, 78/255, 101/255), s=30, label="Control Points",alpha=1)

# 绘制控制多边形
control_points_x = [p[0] for p in control_points]
control_points_y = [p[1] for p in control_points]
# plt.plot(control_points_x, control_points_y, 'r--', label="Control Polygon")  # 使用红色虚线表示控制多边形
plt.plot(control_points_x, control_points_y, color=(201/255, 78/255, 101/255), linestyle='-', linewidth=1.5, label="Control Polygon")

# plt.title("B-Spline Degree=5 in 2D")
plt.xlabel("X",fontsize=20)
plt.ylabel("Y",fontsize=20)
# plt.legend(fontsize=30)
plt.grid(True)
plt.show()

plt.savefig('normal_2D_4degree.pdf')