import numpy as np
import plotly.graph_objects as go
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

# 使用 Plotly 创建 3D 图形
fig = go.Figure()

# 绘制 B-Spline 曲线
fig.add_trace(go.Scatter3d(
    x=y_vals[:, 0], y=y_vals[:, 1], z=y_vals[:, 2], 
    mode='lines', name='B-Spline Degree=4', line=dict(color='blue', width=4)
))

# 绘制控制点
fig.add_trace(go.Scatter3d(
    x=[p[0] for p in control_points], y=[p[1] for p in control_points], z=[p[2] for p in control_points], 
    mode='markers', name='Control Points', marker=dict(color='red', size=6)
))

# 绘制控制多边形
fig.add_trace(go.Scatter3d(
    x=[p[0] for p in control_points], y=[p[1] for p in control_points], z=[p[2] for p in control_points], 
    mode='lines', name='Control Polygon', line=dict(color='red', dash='dash', width=2)
))

# 更新布局
fig.update_layout(
    # title="3D B-Spline Curve with Control Polygon",
    scene=dict(
        xaxis=dict(title='X', range=[min(y_vals[:, 0]) - 0.1, max(y_vals[:, 0]) + 0.1]),
        yaxis=dict(title='Y', range=[min(y_vals[:, 1]) - 0.1, max(y_vals[:, 1]) + 0.1]),
        zaxis=dict(title='Z', range=[min(y_vals[:, 2]) - 0.1, max(y_vals[:, 2]) + 0.1]),
        camera=dict(
            eye=dict(x=1.5, y=1.5, z=1.5)  # 调整视角，确保图形显示完全
        )
    ),
    margin=dict(l=0, r=0, b=0, t=40)  # 增加边距，确保图形不被裁切
)

# 显示图形
fig.show()

# 保存为 PNG
fig.write_image("hard_3D_4degree_plotly.png")
