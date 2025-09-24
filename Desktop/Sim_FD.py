import numpy as np
import math
import matplotlib.pyplot as plt

# ============ 参数与初始值 ============
dt = 1/30

# 鱼初始位置 / 朝向
x_f = [2.0]
y_f = [0.0]
the_f = [0.0]

# 海豚初始位置 / 朝向
x_d = [0.0]
y_d = [0.0]
the_d = [0.0]

# 鱼和海豚的角速度
w_f = [0.0]
w_d = [0.0]

# 鱼、海豚速度
v_f = 9.1549 * 0.4
v_d = 10.1951 * 0.4

# 区域边界
x_min, x_max = -4, 4
y_min, y_max = -4, 4
thre_wall = abs(x_max)/1

def dist_to_wall(x, y, xmin, xmax, ymin, ymax):
    """返回点 (x,y) 到矩形区域 [xmin,xmax]x[ymin,ymax] 边界的最小距离。"""
    dx = min(abs(x - xmin), abs(xmax - x))
    dy = min(abs(y - ymin), abs(ymax - y))
    d_wall = ((x - abs(x_max))**2 + (y-abs(ymax))**2)**0.5
    return min(dx, dy, d_wall)

# ============ 动力学参数举例 ============
# 鱼的转弯动力学
alpha_f  = 14.9719*4*0.25
sigma_f  = 21.4047*1.6*0.5
kappa_f  = 34.7634*20*0.4

# 海豚的转弯动力学
alpha_d = 14.8444*0.4
sigma_d = 13.5023*0.4
kappa_d = 25.7994*0.35
mu      = 2*0.35
kappa_e = 0.0

# fountain event 相关
delta = 0.0001
ds    = 0.5

# 鱼、海豚各自的墙体避让系数 a、b
a_wall_fish    = 2
b_wall_fish    = 4
a_wall_dolphin = 2
b_wall_dolphin = 4

# 仿真步数
num_steps = 1000

# ============ 画图 ============
fig, ax = plt.subplots()
plt.ion()

distances = []

for k in range(num_steps):
    #====================================================#
    #              1) 鱼的位姿更新 (X,Y,Θ)               #
    #====================================================#
    # 根据当前朝向 the_f[k] 计算位移
    x_f_new = x_f[k] + v_f * np.cos(the_f[k]) * dt
    y_f_new = y_f[k] + v_f * np.sin(the_f[k]) * dt
    the_f_new = the_f[k] + w_f[k] * dt
    the_f_new = the_f_new % (2*np.pi)

    # 计算离墙最小距离 df
    df = dist_to_wall(x_f_new, y_f_new, x_min, x_max, y_min, y_max)

    # 如果离墙距离 < 2，则根据公式 (8) 做“指数避让”
    if df < thre_wall:
        phi_f = the_f_new
        avoidance_turn_f = a_wall_fish * np.sign(phi_f) * np.exp(-b_wall_fish * df)
        the_f_new += avoidance_turn_f

    # 更新鱼
    x_f.append(x_f_new)
    y_f.append(y_f_new)
    the_f.append(the_f_new % (2*np.pi))

    #====================================================#
    #            2) 海豚的位姿更新 (X,Y,Θ)              #
    #====================================================#
    x_d_new = x_d[k] + v_d * np.cos(the_d[k]) * dt
    y_d_new = y_d[k] + v_d * np.sin(the_d[k]) * dt
    the_d_new = the_d[k] + w_d[k] * dt
    the_d_new = the_d_new % (2*np.pi)

    # 计算离墙最小距离 dd
    dd = dist_to_wall(x_d_new, y_d_new, x_min, x_max, y_min, y_max)

    # 如果离墙距离 < 2，也用公式 (8) 避让
    if dd < thre_wall:
        phi_d = the_d_new
        avoidance_turn_d = a_wall_dolphin * np.sign(phi_d) * np.exp(-b_wall_dolphin * dd)
        the_d_new += avoidance_turn_d

    # 更新海豚
    x_d.append(x_d_new)
    y_d.append(y_d_new)
    the_d.append(the_d_new % (2*np.pi))

    #====================================================#
    #        3) 计算鱼与海豚的距离，并做fountain事件      #
    #====================================================#
    cur_dis = np.sqrt((x_f_new - x_d_new)**2 + (y_f_new - y_d_new)**2)
    distances.append(cur_dis)

    # 计算鱼对海豚、海豚对鱼的相对角
    dir_f = np.array([np.cos(the_f[k]), np.sin(the_f[k])])
    dir_d = np.array([np.cos(the_d[k]), np.sin(the_d[k])])
    v2 = np.array([x_d[k] - x_f[k], y_d[k] - y_f[k]])
    v3 = np.array([x_f[k] - x_d[k], y_f[k] - y_d[k]])

    Phi_f2d = np.arctan2(np.linalg.det([dir_f, v2]), np.dot(dir_f, v2))
    Phi_f2d = (Phi_f2d + np.pi) % (2 * np.pi) - np.pi

    Phi_d2f = np.arctan2(np.linalg.det([dir_d, v3]), np.dot(dir_d, v3))
    Phi_d2f = (Phi_d2f + np.pi) % (2 * np.pi) - np.pi

    # 计算 fountain 触发系数 c (对鱼)
    c = 1/(1 + np.exp(delta * (cur_dis - ds))) * (1 - np.sign(np.cos(Phi_f2d))) / 2

    # ---- 关键修改：如果鱼离墙太近，就不触发 fountain maneuver ----
    if df < thre_wall:
        c = 0.0

    #====================================================#
    #     4) 反馈控制 + 随机扩散 (更新 w_f, w_d )         #
    #====================================================#
    wf_s = kappa_f * c * Phi_f2d
    wd_s = kappa_d * np.tanh(mu * Phi_d2f) + kappa_e * Phi_d2f

    db1 = np.sqrt(dt)*np.random.randn()
    db2 = np.sqrt(dt)*np.random.randn()

    w_f_new = w_f[k] + (-alpha_f * w_f[k] + wf_s)*dt + sigma_f * db1
    w_d_new = w_d[k] + (-alpha_d * w_d[k] + wd_s)*dt + sigma_d * db2

    w_f.append(w_f_new)
    w_d.append(w_d_new)

    #====================================================#
    #                  5) 绘图                          #
    #====================================================#
    ax.clear()
    ax.plot(x_f, y_f, 'b-', label='Fish')
    ax.plot(x_d, y_d, 'r-', label='Dolphin')
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.legend()
    plt.pause(0.05)

plt.ioff()
plt.show()
