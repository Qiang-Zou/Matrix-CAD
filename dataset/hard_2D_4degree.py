# import numpy as np
# import plotly.graph_objects as go
# from scipy.interpolate import BSpline

# # Degree=5的knots和控制点
# knots = [0, 0, 0, 0, 0, 
#          0.0385, 0.0769, 0.1154, 0.1538, 0.1923, 0.2308, 0.2692, 0.3077, 0.3462, 0.3846,
#          0.4231, 0.4615, 0.5, 0.5385, 0.5769, 0.6154, 0.6538, 0.6923, 0.7308, 0.7692,
#          0.8077, 0.8462, 0.8846, 0.9231, 0.9615,
#          1, 1, 1, 1, 1]  # 尾部节点 (p+1)

# control_points = [
#     [0.0, 0.0],
#     [0.03, 0.3],
#     [0.06, 0.6],
#     [0.1, 1.0],
#     [0.13, 0.6],
#     [0.16, 0.3],
#     [0.2, 0.0],
#     [0.23, 0.3],
#     [0.26, 0.6],
#     [0.3, 1.0],
#     [0.33, 0.6],
#     [0.36, 0.3],
#     [0.4, 0.0],
#     [0.43, 0.3],
#     [0.46, 0.6],
#     [0.5, 1.0],
#     [0.53, 0.6],
#     [0.56, 0.3],
#     [0.6, 0.0],
#     [0.63, 0.3],
#     [0.66, 0.6],
#     [0.7, 1.0],
#     [0.73, 0.6],
#     [0.76, 0.3],
#     [0.8, 1.0],
#     [0.83, 0.6],
#     [0.86, 0.3],
#     [0.9, 0.0],
#     [0.93, 0.3],
#     [0.96, 0.6]
# ]

# # Degree 5 B-Spline曲线
# degree = 4
# spl = BSpline(knots, control_points, degree)

# # 在[0, 1]区间内生成1000个样本点
# x_vals = np.linspace(0, 1, 1000)
# y_vals = spl(x_vals)

# # 使用 Plotly 创建 2D B-Spline 图形
# fig = go.Figure()

# # 添加 B-Spline 曲线
# fig.add_trace(go.Scatter(x=y_vals[:, 0], y=y_vals[:, 1], mode='lines', name='B-Spline Degree=4', line=dict(color='blue', width=2)))

# # 添加控制点
# control_points_x = [p[0] for p in control_points]
# control_points_y = [p[1] for p in control_points]
# fig.add_trace(go.Scatter(x=control_points_x, y=control_points_y, mode='markers', name='Control Points', marker=dict(color='red', size=8)))

# # 添加控制多边形
# fig.add_trace(go.Scatter(x=control_points_x, y=control_points_y, mode='lines', name='Control Polygon', line=dict(color='red', dash='dash')))

# # 更新布局以模拟 MATLAB 风格
# fig.update_layout(
#     # title="B-Spline Degree=4 in 2D",
#     xaxis_title="X",
#     yaxis_title="Y",
#     showlegend=True,
#     template="plotly_white",  # 选择简洁的白色模板
#     plot_bgcolor='white',  # 背景颜色为白色
#     xaxis=dict(showgrid=True, gridcolor='gray', zeroline=True),
#     yaxis=dict(showgrid=True, gridcolor='gray', zeroline=True)
# )

# # 保存图像为 PNG
# fig.write_image("hard_2D_4degree.png")
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
        [0.0, 0.0],
    [0.03, 0.3],
    [0.06, 0.6],
    [0.1, 1.0],
    [0.13, 0.6],
    [0.16, 0.3],
    [0.2, 0.0],
    [0.23, 0.3],
    [0.26, 0.6],
    [0.3, 1.0],
    [0.33, 0.6],
    [0.36, 0.3],
    [0.4, 0.0],
    [0.43, 0.3],
    [0.46, 0.6],
    [0.5, 1.0],
    [0.53, 0.6],
    [0.56, 0.3],
    [0.6, 0.0],
    [0.63, 0.3],
    [0.66, 0.6],
    [0.7, 1.0],
    [0.73, 0.6],
    [0.76, 0.3],
    [0.8, 1.0],
    [0.83, 0.6],
    [0.86, 0.3],
    [0.9, 0.0],
    [0.93, 0.3],
    [0.96, 0.6]
]

print("Length of knots:", len(knots))
print("Shape of control_points:", np.array(control_points).shape)

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

plt.savefig('hard_2D_4degree.pdf')