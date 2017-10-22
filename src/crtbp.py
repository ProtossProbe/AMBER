# -*- coding: UTF-8 -*-

"""
圆型限制性三体问题模块。
"""
import progressbar
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint, ode
from scipy.interpolate import interp1d

mu_m = 0.0121506683
mu_e = 0.0000030359
L1_m = 0.836914718893
L2_m = 1.155682483479
L3_m = -1.00506268026
L1_e = 0.989990937177
L2_e = 1.010070187506
L3_e = -1.00000101196
Radius_e = 6371
Radius_m = 1737.4
LEO = 167
LLO = 150
v_LEO = 7.805
v_LLO = 1.612
v_factor_m = 1.023
v_factor_e = 29.784
t_factor_m = 2.36059488e6/(2*np.pi)
t_factor_e = 3.15581184e7/(2*np.pi)


zhfont = matplotlib.font_manager.FontProperties(
    fname='/Users/huangyukun/Library/Fonts/SourceHanSansCN-Regular.otf')
data_dir = 'data/'


def generate_amatrix(x, y, z, mu=mu_m):
    """
    生成方程的系数矩阵A
    """
    R1squared = y**2 + z**2 + (x + mu)**2
    R2squared = y**2 + z**2 + (-1. + x + mu)**2

    uxx = (1. + 3. * mu * (-1. + x + mu)**2. * R2squared**(-5. / 2.) - mu *
           R2squared**(-3. / 2.) + 3. * (1. - mu) * (x + mu)**2. *
           R1squared**(-5. / 2.) - (1. - mu) * R1squared ** (-3. / 2.))

    uxy = (3. * y * mu * (-1. + x + mu) * R2squared**(-5. / 2.) +
           3. * y * (1. - mu) * (x + mu) * R1squared**(-5. / 2.))

    uxz = (3. * z * mu * (-1. + x + mu) * R2squared**(-5. / 2.) +
           3. * z * (1. - mu) * (x + mu) * R1squared**(-5. / 2.))

    uyy = (1. + 3. * y**2. * mu * R2squared**(-5. / 2.) - mu *
           R2squared**(-3. / 2.) + 3. * y**2. * (1. - mu) *
           R1squared**(-5. / 2.) - (1. - mu) * R1squared**(-3. / 2.))

    uyz = (3. * y * z * mu * R2squared**(-5. / 2.) +
           3. * y * z * (1. - mu) * R1squared**(-5. / 2.))

    uzz = (3. * z**2. * mu * R2squared**(-5. / 2.) - mu *
           R2squared**(-3. / 2.) + 3. * z**2. * (1. - mu) *
           R1squared**(-5. / 2.) - (1. - mu) * R1squared**(-3. / 2.))

    return np.array([[0, 0, 0, 1., 0, 0], [0, 0, 0, 0, 1., 0],
                     [0, 0, 0, 0, 0, 1.], [uxx, uxy, uxz, 0, 2., 0],
                     [uxy, uyy, uyz, -2., 0, 0], [uxz, uyz, uzz, 0, 0, 0]])

def jacobi_inte(y, mu=mu_m,option="normal"):
    """
    计算雅可比积分
    """
    x, y ,z ,dotx, doty, dotz = y[0], y[1], y[2], y[3], y[4], y[5]
    U = 0.5*(x**2 + y**2) + 0.5*(1 - mu)*mu + mu/np.sqrt(y**2 + z**2 + (-1 + \
        x + mu)**2) + (1 - mu)/np.sqrt(y**2 + z**2 + (x + mu)**2)
    if option is "normal":
        return 2*U-dotx**2-doty**2-dotz**2
    else:
        return 2*U-dotx**2-doty**2-dotz**2-mu*(1-mu)

def mani_jacobi(manifolds, mu=mu_m):

    n = len(manifolds)
    j1=0
    for traj in manifolds:
        temp = np.mean(jacobi_inte(traj[:,1:].T, mu))
        j1 = j1 + temp
    return j1 / n

def plot_jacobi_range(C, delta, mu=mu_m):
    """
    画出雅可比积分的等高线图
    """
    # 图绘制的范围
    range = 1.5

    # 图绘制的中心（x轴）
    center = 0

    # 计算各个坐标下的雅各比积分
    x = np.arange(-range + center, range + center, delta)
    y = np.arange(-range, range, delta)
    X, Y = np.meshgrid(l)
    Z = jacobi_inte([X, Y, 0, 0, 0, 0], mu)

    # 绘出等高线图
    fig, ax = plt.subplots()
    p = ax.contour(Z, levels = C, extent=[x[0], x[-1],  \
        y[0], y[-1]])
    ax.plot([1-mu], [0], 'oy')
    ax.plot([-mu], [0], 'ob')
    an1 = ax.annotate(u'月球', xy=(1-mu,0.05), ha='center',fontproperties=zhfont)
    an2 = ax.annotate(u'地球', xy=(-mu,0.05), ha='center',fontproperties=zhfont)
    an3 = ax.annotate(u'禁行域', xy=(0,1), ha='center',fontproperties=zhfont)
    an4 = ax.annotate(u'希尔域', xy=(0,0.4), fontsize=8,ha='center', fontproperties=zhfont)
#    ax.plot([L1_m,L2_m,L3_m], [0,0,0], 'xk')
    #cb = fig.colorbar(p, ax=ax)

    # 额外参数
    ax.set_aspect('equal')
    #cb.set_label(u'雅可比积分', fontproperties=zhfont)
    ax.set_xlabel(u'x轴', fontproperties=zhfont)
    ax.set_ylabel(u'y轴', fontproperties=zhfont)
    plt.savefig('zero_curve.pdf')

def crtbpfunc(y, t, mu=mu_m):
    """
    圆型限制性三体问题的方程，在六维状态空间中给出，适用于odeint
    """

    y1, y2, y3, y4, y5, y6 = y[0], y[1], y[2], y[3], y[4], y[5]

    dy1 = y4
    dy2 = y5
    dy3 = y6
    dy4 = -(1-mu)*(y1+mu)/((y1+mu)**2+y2**2+y3**2)**(3./2.)-mu*(y1-1+mu)/((y1 \
            -1+mu)**2+y2**2+y3**2)**(3./2.)+y1+2*y5
    dy5 = -(1-mu)*y2/((y1+mu)**2+y2**2+y3**2)**(3./2.)-mu*y2/((y1-1+mu)**2+y2 \
            **2+y3**2)**(3./2.)+y2-2*y4
    dy6 = -((y3*mu)/(y2**2+y3**2+(-1+y1+mu)**2)**(3./2.))-(y3*(1-mu))/(y2**2+ \
            y3**2+(y1+mu)**2)**(3./2.)

    return [dy1, dy2, dy3, dy4, dy5, dy6]


def crtbpfunc_lin(y, t, mu=mu_m):
    """
    圆型限制性三体问题的线性化方程，在六维状态空间中给出，适用于odeint
    """
    A = generate_amatrix(L2_m, 0, 0, mu)
    dy = A.dot(y)
    return dy

def crtbpfunc_tran(t, y, mu=mu_m):
    """
    圆型限制性三体问题的方程，在六维状态空间中给出，适用于ode
    """

    y1, y2, y3, y4, y5, y6 = y[0], y[1], y[2], y[3], y[4], y[5]

    dy1 = y4
    dy2 = y5
    dy3 = y6
    dy4 = -(1-mu)*(y1+mu)/((y1+mu)**2+y2**2+y3**2)**(3./2.)-mu*(y1-1+mu)/((y1 \
            -1+mu)**2+y2**2+y3**2)**(3./2.)+y1+2*y5
    dy5 = -(1-mu)*y2/((y1+mu)**2+y2**2+y3**2)**(3./2.)-mu*y2/((y1-1+mu)**2+y2 \
            **2+y3**2)**(3./2.)+y2-2*y4
    dy6 = -((y3*mu)/(y2**2+y3**2+(-1+y1+mu)**2)**(3./2.))-(y3*(1-mu))/(y2**2+ \
            y3**2+(y1+mu)**2)**(3./2.)

    return [dy1, dy2, dy3, dy4, dy5, dy6]

def a_transform(A):
    """
    将矩阵A转化为36阶方阵用于积分
    """
    At = np.array([]).reshape(0,36)
    for i in np.arange(6):
        At1 = np.array([]).reshape(6,0)
        for j in np.arange(6):
            temp = np.eye(6) * A[i, j]
            At1 = np.hstack((At1, temp))
        At = np.vstack((At,At1))
    return At

def crtbpfunc42(y, t, mu=mu_m):
    """
    用于微分矫正法生成周期轨道的42个微分方程，适用于odeint
    """
    dy = crtbpfunc(y[0:6], t, mu)
    phi = y[6:42]
    At = a_transform(generate_amatrix(y[0], y[1], y[2], mu))
    dphi = np.dot(At, phi)
    return np.append(dy, dphi)

def crtbpfunc42_tran(t, y, mu=mu_m):
    """
    用于微分矫正法生成周期轨道的42个微分方程，适用于ode
    """
    dy = crtbpfunc(y[0:6], t, mu)
    phi = y[6:42]
    At = a_transform(generate_amatrix(y[0], y[1], y[2], mu))
    dphi = np.dot(At, phi)
    return np.append(dy, dphi)

def jacobi_error(y, mu=mu_m, plot=False):
    """
    给定一系列的轨道状态，分别计算雅可比积分，返回雅可比积分的波动范围以及平均雅可比积分
    """
    # 计算雅可比积分
    jacobi = jacobi_inte(y.T, mu)
    Meanjacobi = np.mean(jacobi)

    # 计算雅可比积分的误差
    error = jacobi - Meanjacobi
    Maxerror = np.max(error)
    Minerror = np.min(error)

    # 在给定条件下画图
    if plot:
        fig, axes = plt.subplots(figsize=(4,4))
        axes.plot(np.arange(len(error)), error, 'r')
        axes.set_xlabel('index')
        axes.set_ylabel('error')
        axes.set_title('error fluctuation')

    return np.array([Maxerror, Minerror, Meanjacobi])

def orbit_plot(y, range=0.5, target='moon', option='2d'):
    """
    给定一系列的轨道状态，画出三维空间下的轨道图
    """
    if option is '3d':
        # 初始化3D绘图模块
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # 指定绘图范围
        ax.set_zlim3d(-range,range)
        ax.set_xlim3d(1-range, 1+range)
        ax.set_ylim3d(-range, range)

        # 画出轨道
        if type(y) is np.ndarray:
            ax.plot(y[:,0], y[:,1], y[:,2], 'g',alpha=0.7)
        elif type(y) is list:
            for i in y:
                ax.plot(i[:,1], i[:,2], i[:,3], 'g',alpha=0.7)

        # 两个天体与拉格朗日点
        if target is 'moon':
            plot_primaries(ax, option='em3d')
        elif target is 'earth':
            plot_primaries(ax, option='se3d')

    elif option is '2d':
        # 初始化2D绘图模块
        fig,ax = plt.subplots()

        # 指定绘图范围
        ax.set_xlim(1-range, 1+range)
        ax.set_ylim(-range, range)

        # 画出轨道
        if type(y) is np.ndarray:
            ax.plot(y[:,0], y[:,1],alpha=0.7, linewidth=1.5)
        elif type(y) is list:
            for i in y:
                ax.plot(i[:,1], i[:,2],alpha=0.7, linewidth=1.5)

        # 两个天体与拉格朗日点
        if target is 'moon':
            LEO_orbit = plt.Circle((-mu_m,0),0.017008,color='g',fill=False,
                                   ls='dashed')
            LLO_orbit = plt.Circle((1-mu_m,0),0.00491,color='y',fill=False,
                                   ls='dashed')
            fig.gca().add_artist(LEO_orbit)
            fig.gca().add_artist(LLO_orbit)
            plot_primaries(ax, option='em2d')
        elif target is 'earth':
            LEO_orbit = plt.Circle((1-mu_e,0),4.3703e-05,color='g',fill=False,
                                   ls='dashed')
            fig.gca().add_artist(LEO_orbit)
            plot_primaries(ax, option='se2d')

    # 图片展示
    #plt.savefig('L2_unstable1.pdf')

def save_orbit_data(filename, y0):
    """
    存储新的轨道初值信息到.npy文件中
    """
    y_data = np.load(filename)
    y_data = np.vstack((y_data,y0))
    np.save(filename, y_data)

def odeint42(y0, dt, acc=1e-12, target='moon', option='half',t_end=10):
    """
    对给定的初始条件进行积分，返回状态、时间以及Phi矩阵
    """
    # 初始化Phi矩阵
    phi0 = np.eye(6).reshape(1,36)

    # 积分初始条件设置
    init = np.append(y0, phi0)
    t0 = 0
    t1 = t_end
    value_t = np.array([0.])
    value_y = init
    temp = init

    # 选择积分方法
    backend = "dopri5"
    #backend = "dop853"
    #backend = "lsoda"
    #backend = "vode"

    # 设置积分器
    r = ode(crtbpfunc42_tran).set_integrator(backend)
    if target is 'moon':
        r.set_initial_value(init, t0)
    elif target is 'earth':
        r.set_initial_value(init, t0).set_f_params(mu_e)

    # 如果选项是 half，则积分一直到穿过x-z平面时，靠近积分终点时逐步缩小步长
    if option is 'half':
        while r.successful() and r.t < t1:
            r.integrate(r.t+dt, step=True)

            # 如果检测到未穿过x-z平面，则储存积分信息
            if temp[1] * r.y[1] >= 0:
                value_t = np.append(value_t, r.t)
                value_y = np.vstack((value_y, r.y.T))
                temp = r.y

            # 如果穿过x-z平面，并且精度不符合要求，则缩小步长，并退回到积分的上一步
            elif np.abs(r.y.T[1]) > acc:
                r.set_initial_value(temp, r.t-dt)
                dt = dt / 10.0
            # 如果穿过x-z平面，并且精度符合要求，则积分结束
            elif np.abs(r.y.T[1]) < acc:
                break

    # 如果选项是 1cycle，则积分一个周期，步长不变
    elif option is '1cycle':
        flag = 0

        # 积分直到一个周期
        while r.successful() and r.t < t1 and flag < 2:
            r.integrate(r.t+dt, step=True)
            value_t = np.append(value_t, r.t)
            value_y = np.vstack((value_y, r.y.T))
            if temp[1] * r.y[1] < 0:
                flag = flag + 1
            temp = r.y

    elif option is 'time':
        while r.successful() and np.abs(r.t - t1) > 1e-12:
            if r.t - t1 < 0:
                r.integrate(r.t+dt, step=True)
                value_t = np.append(value_t, r.t)
                value_y = np.vstack((value_y, r.y.T))
                temp = r.y
            elif r.t - t1 > 0:
                r.set_initial_value(temp, r.t-dt)
                dt = dt / 10.0

    # 返回积分结果
    y = value_y[:,0:6]
    phi = value_y[:,6:42]
    phi = phi.reshape(phi.shape[0],6,6)
    phi_end = phi[-1,:]

    return y, value_t, phi, phi_end

def cal_acc(y,model='crtbp',mu=mu_m):
    """
    给定六维相空间中的一点，计算该点的加速度大小
    """
    y1, y2, y3, y4, y5, y6 = y[0], y[1], y[2], y[3], y[4], y[5]
    if model is 'crtbp':
        xdoubledot = -(1-mu)*(y1+mu)/((y1+mu)**2+y2**2+y3**2)**(1.5)-mu*(y1-1+mu)/ \
                    ((y1-1+mu)**2+y2**2+ y3**2)**(1.5)+y1+2*y5
        ydoubledot = -(1-mu)*y2/((y1+mu)**2+y2**2+y3**2)**(1.5)-mu*y2/((y1-1+mu)**2+y2 \
                    **2+y3**2)**(1.5)+y2-2*y4
        zdoubledot = -((y3*mu)/(y2**2+y3**2+(-1+y1+mu)**2)**(1.5))-(y3*(1-mu))/(y2 \
                    **2+y3**2+(y1+mu)**2)**(1.5)
    return np.array([xdoubledot,ydoubledot,zdoubledot]).T

def diff_correct(y0, y, phi, option="fix_x", mu=mu_m):
    """
    微分矫正法返回修正值
    """
    # 初始化变量
    y1, y2, y3, y4, y5, y6 = y[0], y[1], y[2], y[3], y[4], y[5]
    y0_1, y0_2, y0_3, y0_4, y0_5, y0_6 = y0[0], y0[1], y0[2], y0[3], y0[4], y0 \
    [5]

    # 计算xdoubledot和zdoubledot
    temp = cal_acc(y,mu=mu)
    xdoubledot, zdoubledot = temp[0],temp[2]

    # 计算雅可比积分的偏导数
    Cy0dot = -2*y0_5
    Cx0 = 2*(-(1-mu)*(y0_1+mu)/((y0_1+mu)**2+y0_2**2+y0_3**2)**(1.5)-mu* \
         (y0_1-1+mu)/((y0_1-1+mu)**2+y0_2**2+y0_3**2)**(1.5)+y0_1)
    Cz0 = 2*(-((y0_3*mu)/(y0_2**2+y0_3**2+(-1+y0_1+mu)**2)**(1.5))-(y0_3* \
         (1-mu))/(y0_2**2+y0_3**2+(y0_1+mu)**2)**(1.5))

    xy = xdoubledot/y5
    zy = zdoubledot/y5

    A = np.eye(6)
    A[1,1] = -1
    A[3,3] = -1
    A[5,5] = -1

#    B1 = np.array([[phi[3,2],phi[3,4]],[phi[5,2],phi[5,4]]])-1./y5*np.dot \
#        (np.array([[xdoubledot,zdoubledot]]).T,np.array([[phi[1,2],phi[1,4]]]))
#    B2 = np.array([[phi[3,0],phi[3,4]],[phi[5,0],phi[5,4]]])-1./y5*np.dot \
#        (np.array([[xdoubledot,zdoubledot]]).T,np.array([[phi[1,0],phi[1,4]]]))

    P = np.zeros([3,3])
    P[0,0] = phi[3,0] - phi[1,0] * xy
    P[0,1] = phi[3,2] - phi[1,2] * xy
    P[0,2] = phi[3,4] - phi[1,4] * xy
    P[1,0] = phi[5,0] - phi[1,0] * zy
    P[1,1] = phi[5,2] - phi[1,2] * zy
    P[1,2] = phi[5,4] - phi[1,4] * zy
    P[2,0] = Cx0
    P[2,1] = Cz0
    P[2,2] = Cy0dot

#    M1 = np.array([[phi[0,0],phi[0,2],phi[0,3],phi[0,5]],[phi[2,0],phi[2,2], \
#        phi[2,3],phi[2,5]],[phi[3,0],phi[3,2],phi[3,3],phi[3,5]],[phi[5,0], \
#        phi[5,2],phi[5,3],phi[5,5]]])
#    M2 = 1/Cy0dot * np.dot(np.array([[phi[0,4],phi[2,4],phi[3,4],phi[5,4]]]). \
#        reshape(4,1),np.array([[Cx0,Cz0,0,0]]))
#    M4 = np.array([[phi[1,0],phi[1,2],phi[1,3],phi[1,5]]]) - phi[1,4]/Cy0dot * \
#        np.array([[Cx0,Cz0,0,0]])
#    M3 = 1/y5 * np.dot(np.array([[0,0,xdoubledot,zdoubledot]]).reshape(4,1),M4)
#
#    Monodromy1 = M1 - M2 - M3

    # 计算单值矩阵M
    mono = A.dot(np.linalg.inv(phi)).dot(A).dot(phi)

    # 在不同方法下返回不同的修正值
    if option == "fix_x":
        result_fix_x = np.linalg.solve(P[0:2,1:3],np.array([-y4,-y6]))
        output = np.array([0, 0, result_fix_x[0], 0, result_fix_x[1], 0])
    elif option == "fix_z":
        result_fix_z = np.linalg.solve(P[0:2,0::2],np.array([-y4,-y6]))
        output = np.array([result_fix_z[0], 0, 0, 0, result_fix_z[1], 0])
    elif option == "fix_c":
        result_fix_c = np.linalg.solve(P,np.array([-y4,-y6,0]))
        output = np.array([result_fix_c[0], 0, result_fix_c[1], 0, result_fix_c[2], 0])

    # 计算误差值
    dotx_error = y4
    dotz_error = y6
    error = np.mean(np.abs(dotx_error)+np.abs(dotz_error))

    return output, error, mono

def multi_shooting_step_1(y_p, y, phi, option="normal", mu=mu_m):
    """
    多步打靶法的第一步：固定每一个分段的位置，修正速度

    公式参考文献：Numerical determination of Lissajous trajectories in the
    restricted three-body problem
    """
    # 初始化变量
    y1, y2, y3, y4, y5, y6 = y[0], y[1], y[2], y[3], y[4], y[5]
    y_p1, y_p2, y_p3, y_p4, y_p5, y_p6 = y_p[0], y_p[1], y_p[2], y_p[3], \
                                        y_p[4], y_p[5]

    L = phi[0:3,3:]
    L = np.hstack((L,np.array([[y4],[y5],[y6]])))
    deltax_p = y_p1 - y1
    deltay_p = y_p2 - y2
    deltaz_p = y_p3 - y3
    b = np.array([[deltax_p,deltay_p,deltaz_p]]).T
    L_temp = np.linalg.inv(L.dot(L.T))
    u = L.T.dot(L_temp).dot(b).T

    y_corr = np.hstack((np.zeros((1,3)),u[:,0:3]))
    time_corr = u[:,-1]

    return y_corr, time_corr

def phi_slice(phi):
    A = phi[0:3,0:3]
    B = phi[0:3,3:]
    C = phi[3:,0:3]
    D = phi[3:,3:]
    return A, B, C, D

def mmatrix(y_minus,y_plus,phi_0,phi_1,mu=mu_m):
    phi_1 = np.linalg.inv(phi_1)
    A0,B0,C0,D0 = phi_slice(phi_0)
    Af,Bf,Cf,Df = phi_slice(phi_1)
    a_minus = cal_acc(y_minus).reshape(3,1)
    a_plus = cal_acc(y_plus).reshape(3,1)
    v_minus = y_minus[3:].reshape(3,1)
    v_plus = y_plus[3:].reshape(3,1)

    M0 = D0.dot(np.linalg.inv(B0)).dot(A0) - C0
    Mt0 = a_minus - D0.dot(np.linalg.inv(B0)).dot(v_minus)
    Mp = Df.dot(np.linalg.inv(Bf)) - D0.dot(np.linalg.inv(B0))
    Mtp = D0.dot(np.linalg.inv(B0)).dot(v_minus) - Df.dot(np.linalg.inv(Bf)) \
            .dot(v_plus) + a_plus - a_minus
    Mf = Cf - Df.dot(np.linalg.inv(Bf)).dot(Af)
    Mtf = Df.dot(np.linalg.inv(Bf)).dot(v_plus) - a_plus

    M = np.hstack((M0,Mt0,Mp,Mtp,Mf,Mtf))
    return M

def generate_mmatrix(y_aft_minus,y_aft_plus,Phi,mu=mu_m):
    n = Phi.shape[0]
    M = np.array([]).reshape(0,4*n+4)
    Phi_0 = Phi[0:-1]
    Phi_1 = Phi[1:]
    index = 0
    for y_minus,y_plus,phi_0,phi_1 in zip(y_aft_minus,y_aft_plus,Phi_0,Phi_1):
        m = mmatrix(y_minus,y_plus,phi_0,phi_1,mu)
        m = np.hstack((np.zeros((3,index*4)),m,np.zeros((3,4*n+4-index*4-12))))
        M = np.vstack((M,m))
        index = index + 1
    return M

def multi_shooting_step_2(y_aft_minus,y_aft_plus,Phi,mu=mu_m,):
    n = Phi.shape[0]
    M = generate_mmatrix(y_aft_minus,y_aft_plus,Phi,mu=mu_m)
    delta_y_aft = y_aft_plus-y_aft_minus
    delta_v = delta_y_aft[:,3:].reshape(3*(n-1),1)
    v_error = 0
    for v in delta_y_aft[:,3:]:
        v = np.linalg.norm(v)
        v_error = v_error + v
    temp = np.linalg.inv(M.dot(M.T))
    delta_r = -M.T.dot(temp).dot(delta_v)
    delta_r = delta_r.reshape(n+1,4)
    return delta_r,v_error


def cal_eigenvector(mono):
    """
    使用单值矩阵返回初始状态的不稳定特征向量和稳定特征向量
    """
    # 计算单值矩阵的特征值和特征向量
    w,v = np.linalg.eig(mono)

    # 不稳定特征向量是特征值最大的那一个
    unstable_index = np.argmax(w.real)
    v_unstable = v[:,unstable_index]

    # 稳定特征向量是特征值最小的那一个
    stable_index = np.argmin(w.real)
    v_stable = v[:,stable_index]

    return np.array([v_unstable.real, v_stable.real])

def manifolds_init(y, phi, v, Period, target='moon', numberofstep=100):
    """
    使用单值矩阵和周期轨道的初始条件，生成周期轨道不变流形的初始条件
    """
    # 周期轨道上一共有n个数据
    n = y.shape[0]

    # 数据的时间间隔
    dt = Period / n

    # 一个等分间隔的时间
    step = Period/numberofstep/dt

    # 均分的指标
    index = [int(np.floor(i*step)) for i in np.arange(numberofstep)]

    # 得到每个位置的稳定与不稳定向量，以及初始条件
    v_unstable = np.array([phi[i].dot(v[0]) for i in index])
    v_stable = np.array([phi[i].dot(v[1]) for i in index])
    y_new = y[index]

    # 初始化变量
    y_stable1 = y_stable2 = y_unstable1 = y_unstable2 = np.array([]).reshape(0,6)

    # 生成稳定与不稳定流形的初始条件（共4组）
    for i in np.arange(numberofstep):

        # 设置扰动量的大小
        normv = np.linalg.norm(y_new[i,3:6])
        if target is 'moon':
            epsilon = np.abs(100./384400./normv)
        elif target is 'earth':
            epsilon = np.abs(500./1.496e8/normv)
        # 生成稳定与不稳定流形的初始条件
        y_stable1 = np.vstack((y_stable1, y_new[i] + epsilon * v_stable[i]/np \
                .linalg.norm(v_stable[i])))
        y_stable2 = np.vstack((y_stable2, y_new[i] - epsilon * v_stable[i]/np. \
                linalg.norm(v_stable[i])))
        y_unstable1 = np.vstack((y_unstable1, y_new[i] + epsilon * v_unstable \
                [i]/np.linalg.norm(v_unstable[i])))
        y_unstable2 = np.vstack((y_unstable2, y_new[i] - epsilon * v_unstable \
                [i]/np.linalg.norm(v_unstable[i])))

    return y_stable1, y_stable2, y_unstable1, y_unstable2

def manifolds_gene(y_series, option, time=2, acc=1e-11, end='secondary', \
                    mu=mu_m, target='moon', crosstime=1, angle=np.pi/2):
    """
    使用流形的初始条件，生成一组不变流形
    """

    # 设置积分步长
    step = 0.001

    # 通过option选择积分稳定流形还是不稳定流形
    if option is 'u':
        step = step
    elif option is 's':
        step = -step

    # 共n组初始条件
    n = y_series.shape[0]

    # 条件初始化
    t0 = 0
    t1 = time
    manifolds = []
    if end is 'secondary':
        cut = 1-mu
    elif end is 'primary':
        cut = -mu

    # 选择积分方法
    #backend = "dopri5"
    backend = "dop853"
    #backend = "lsoda"
    #backend = "vode"

    # 设置积分器
    r = ode(crtbpfunc_tran).set_integrator(backend, rtol = 1e-9, atol = 1e-11 \
        , nsteps = 10000)

    bar = progressbar.ProgressBar()
    for i in bar(np.arange(n)):
        flag = 1
        # 设置积分初始条件
        value_t = np.array([0.]).reshape(1,1)
        value_y = y_series[i]
        temp = value_y
        if target is 'moon':
            r.set_initial_value(value_y, t0)
        elif target is 'earth':
            r.set_initial_value(value_y, t0).set_f_params(mu_e)
        dt = step
        while r.successful() and np.abs(r.t) < t1:
            r.integrate(r.t+dt, step=True)
            ifcross, error = cross_section(temp, r.y, cut, angle,option='section')
            # 如果没有穿过所需截面，或者穿过次数不符合要求，就正常积分
            if (not ifcross) or flag < crosstime:
                # 如果检测到穿过截面，就增加flag值
                if ifcross:
                    flag = flag + 1
                value_t = np.vstack((value_t, r.t))
                value_y = np.vstack((value_y, r.y.T))
                temp = r.y

            # 如果穿过所需截面，且精度未达到要求，则缩小步长，返回上一步积分
            elif error > acc:
                r.set_initial_value(temp, r.t-dt)
                dt = dt / 10.0

            # 如果穿过所需截面，且精度达到要求，则积分完成
            elif error < acc:
                break

        # 一个初始条件积分完毕后，将其储存到一个数组中
        manifolds.append(np.hstack((value_t,value_y)))
    return manifolds

def cross_section(y_now, y_next, cut, theta=np.pi/2, option='section'):
    """
    通过前后两个状态向量，判断是否已经穿过截面
    """
    pos = cut + 0j
    xy_now = y_now[0] + y_now[1]*1j
    xy_next = y_next[0] + y_next[1]*1j
    angle_now = np.angle(xy_now-pos)
    angle_next = np.angle(xy_next-pos)

    if option is 'section':
        # 如果空间上穿过了那一定穿过了
        if (y_now[0]-cut) * (y_next[0]-cut) < 0:
            error = y_next[0]-cut
            return True, np.abs(error)
        # 如果在接近截面处的地方速度反向了，那也一定穿过了
        elif (y_now[3] * y_next[3]) < 0 and np.abs(y_now[0]-cut)+np.abs(y_next[0] \
            -cut) < 0.005:
            error = (angle_next-theta)*np.linalg.norm(xy_next)
            return True, np.abs(error)

        # 其他情况都是没穿过
        else:
            return False, 0

    elif option is 'angle':
        # 如果空间上穿过了那一定穿过了
        if (angle_now-theta) * (angle_next-theta) < 0:
            error = (angle_next-theta)*np.linalg.norm(xy_next)
            return True, np.abs(error)
        # 如果在接近截面处的地方速度反向了，那也一定穿过了
        elif (y_now[3] * y_next[3]) < 0 and np.abs(y_now[0]-cut)+np.abs(y_next[0] \
            -cut) < 0.005:
            error = (angle_next-theta)*np.linalg.norm(xy_next)
            return True, np.abs(error)

        # 其他情况都是没穿过
        else:
            return False, 0

def poincare_section(manifolds, x_cut=0, plot=True):
    """
    画出与庞加莱截面的交线
    """
    # 一共有多少组流形
    n = len(manifolds)
    y_section = np.array([[]]).reshape(0,6)

    # 最后一个数据当然就是穿过截面时的数据
    for i in np.arange(n):
        y_section = np.vstack((y_section, manifolds[i][-1,1:]))

#    while np.abs(y_section[:,3]).max() > 30.:
#        y_section = np.delete(y_section, (np.abs(y_section[:,3]).argmax()), axis=0)

    # 画图
    if plot:
        fig, ax = plt.subplots()
    #    ax.set_xlim(-0.1, 0.15)
    #    ax.set_ylim(-5, 5)
        ax.plot(y_section[:,1], y_section[:,4], 'r.')
        plt.show()
    return y_section

def orbit_intepolate(query, orbit = 'L1_p', option = 'find_y', mu=mu_m):
    """
    使用插值的方法，根据需求返回对应的量
    """
    # 存档后的文件是这些
    name_group = ['orbit_init_data_L2_halo.npy', 'orbit_init_data_L1_halo.npy', \
            'orbit_init_data_8_shape.npy']

    # 通过orbit选择合适的一簇周期轨道的信息
    if orbit is 'L1_p':
        filename = data_dir + name_group[1]
    elif orbit is 'L2_p':
        filename = data_dir + name_group[0]

    # 读取轨道信息
    y_data = np.load(filename)
    y_data = y_data[y_data[:,0].argsort()]

    # 计算这簇轨道的雅各比积分
    jacobi_data = jacobi_inte(y_data.T, mu)

    # 根据option的不同返回不同的值
    method = 'cubic'
    if option is 'find_y':
        spl = interp1d(y_data[:,0], y_data[:,4], kind=method)
    elif option is 'find_jacobi':
        spl = interp1d(y_data[:,0], jacobi_data, kind=method)
    elif option is 'find_x_with_jacobi':
        spl = interp1d(jacobi_data, y_data[:,0], kind=method)
    elif option is 'find_x_with_y':
        spl = interp1d(y_data[:,4], y_data[:,0], kind=method)

    # 返回结果
    result = spl(query)

#    fig, axes = plt.subplots()
#    axes.plot(y_data[:,0], y_data[:,4], '.r')
#    xs = np.linspace(y_data[:,0].min(), y_data[:,0].max(), 100)
#    axes.plot(xs, spl(xs), 'b')
#    plt.show()
    return result

def find_init_with_same_jacobi(jacobi, mu=mu_m):
    """
    指定雅可比积分，返回符合这个雅各比积分数值的周期轨道初始条件
    """
    x_L1 = orbit_intepolate(jacobi, orbit = 'L1_p', option = 'find_x_with_jacobi')
    x_L2 = orbit_intepolate(jacobi, orbit = 'L2_p', option = 'find_x_with_jacobi')
    y_L1 = orbit_intepolate(x_L1, orbit = 'L1_p', option = 'find_y')
    y_L2 = orbit_intepolate(x_L2, orbit = 'L2_p', option = 'find_y')
    return np.array([x_L1,0,0,0,y_L1,0]), np.array([x_L2,0,0,0,y_L2,0])

def plot_poincare_cross(y_section_u, y_section_s, plot_range, option='2d'):
    """
    使用两个流形的庞加莱截面交线数据，画出曲线的交点，并且标出序号
    """
    # 一个流形有多少个点
    n = y_section_u.shape[0]
    m = y_section_s.shape[0]

    if option is '2d':
        # 画图初始化
        fig, ax = plt.subplots()
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        # 设置画图范围
        ax.set_xlim(plot_range[0], plot_range[1])
        ax.set_ylim(plot_range[2], plot_range[3])

        # 画出两个流形及其序号
        ax.plot(y_section_u[:,1], y_section_u[:,4], 'r.')
        for i,j,k in zip(y_section_u[0:-1:10,1], y_section_u[0:-1:10,4], np. \
            arange(n)):
            ax.annotate(str(k), xy=(i,j), fontsize=8, alpha=0.6)
        ax.plot(y_section_s[:,1], y_section_s[:,4], 'g.')
        for i,j,k in zip(y_section_s[0:-1:10,1], y_section_s[0:-1:10,4], np. \
            arange(m)):
            ax.annotate(str(k), xy=(i,j), fontsize=8, alpha=0.6)
        ax.set_xlabel(r'\textbm{y}')
        ax.set_ylabel(r'\dot{y}')
    elif option is '3d':
        # 画图初始化
        fig = plt.figure()
        ax = Axes3D(fig)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        # 设置画图范围
        ax.set_xlim3d(plot_range[0], plot_range[1])
        ax.set_ylim3d(plot_range[2], plot_range[3])
        ax.set_zlim3d(plot_range[4], plot_range[5])

        # 画出两个流形及其序号
        ax.plot(y_section_u[:,1], y_section_u[:,3], y_section_u[:,4], 'r.')
        ax.plot(y_section_s[:,1], y_section_s[:,3], y_section_s[:,4], 'g.')
        ax.set_xlabel(r'\textbm{y}')
        ax.set_ylabel(r'\dot{x}')
        ax.set_zlabel(r'\dot{y}')
    #plt.show()
    # ax.set_xlabel(r'\textbm{y}')
    # ax.set_ylabel(r'\dot{y}')
    # plt.savefig('intersection.pdf')

def jacobi_xdot(y, jacobi, option='positive', mu=mu_m):
    """
    给定雅可比积分数值和初始条件，求出xdot
    """
    x, y ,z ,doty, dotz = y[0], y[1], y[2], y[3], y[4]

    U = 0.5*(x**2 + y**2) + 0.5*(1 - mu)*mu + mu/np.sqrt(y**2 + z**2 + (-1 + \
        x + mu)**2) + (1 - mu)/np.sqrt(y**2 + z**2 + (x + mu)**2)
    dotx = np.sqrt(-jacobi+2*U-mu*(1-mu)-dotz**2-doty**2)
    if option is 'positive':
        return dotx
    elif option is 'negative':
        return -dotx

def section_delete(y_section, option='moon'):
    """
    将进入天体半径范围的截面条件删除
    """
    moon_radius = 0.00452
    earth_radius = 0.00004264
    if option is 'moon':
        y_section[:][np.abs(y_section[:,1])<moon_radius] = 0
    elif option is 'earth':
        y_section[:][np.abs(y_section[:,1])<earth_radius] = 0
    return y_section

def unit_trans(l, v, t, target='moon', option='dimension'):
    """
    无量纲单位转换到真实空间的单位
    """
    if target is 'moon':
        L_factor = 3.844e5
        T_factor = 2.36059488e6/(2*np.pi)
        V_factor = 1.023
    elif target is 'earth':
        L_factor = 1.496e8
        T_factor = 3.15581184e7/(2*np.pi)
        V_factor = 29.784
    if option is 'dimension':
        l = l * L_factor
        t = t * T_factor
        v = v * V_factor
    elif option is 'nondimension':
        l = l / L_factor
        t = t / T_factor
        v = v / V_factor
    return l, v, t

def R_matrix(t, theta0=0.):
    """
    在不同时刻下的旋转矩阵
    """
    theta_t = t + theta0
    c = np.cos(theta_t)
    s = np.sin(theta_t)
    R11 = np.array([[c,-s,0],[s,c,0],[0,0,1]])
    R12 = np.zeros([3,3])
    R21 = np.array([[-s,-c,0],[c,-s,0],[0,0,0]])
    R22 = np.array([[c,-s,0],[s,c,0],[0,0,1]])
    R = np.vstack((np.hstack((R11,R12)),np.hstack((R21,R22))))
    return R

def frame_trans(y_temp, t, fix_p=-mu_m, theta=0.,option='rota2init'):
    """
    不同坐标系的转换
    """
    y = np.copy(y_temp)
    d = np.array([[fix_p,0,0,0,0,0]]).reshape(6,1)
    n = y.shape[0]
    if option is 'rota2init':
        for i in np.arange(n):
            R = R_matrix(t[i], theta)
            y[i] = R.dot(y[i].reshape(6,1)-d).reshape(6,)
    elif option is 'init2rota':
        for i in np.arange(n):
            R = np.linalg.inv(R_matrix(t[i], theta))
            temp = R.dot(y[i].reshape(6,1)) + d
            y[i] = temp.reshape(6,)
    return y

def primary_trans(y_temp, t_temp, option='m2e'):
    """
    转换主天体
    """
    y = np.copy(y_temp)
    t = np.copy(t_temp)
    L_a = 3.844e5
    T_a = 2.36059488e6/(2*np.pi)
    L_b = 1.496e8
    T_b = 3.15581184e7/(2*np.pi)
    if option is 'm2e':
        y[:,0:3] = y[:,0:3] * L_a / L_b
        y[:,3:] = y[:,3:] * L_a * T_b / (L_b * T_a)
        t = t * T_a / T_b
    elif option is 'e2m':
        y[:,0:3] = y[:,0:3] * L_b / L_a
        y[:,3:] = y[:,3:] * L_b * T_a / (L_a * T_b)
        t = t * T_b / T_a
    return y, t

def em2se(y_temp, t_temp, theta0=0.):
    """
    从地-月系转移到日-地月系
    """
    y, t = np.copy(y_temp), np.copy(t_temp)
    y_tran = frame_trans(y, t, fix_p = -mu_m,option='rota2init')
    y_e,t_e = primary_trans(y_tran, t)
    y_e_tran = frame_trans(y_e, t_e, fix_p = 1-mu_e, theta=theta0, \
                option='init2rota')
    return y_e_tran, t_e

def se2em(y_temp, t_temp, theta0=0.):
    """
    从日-地/月系转移到地-月系
    """
    y, t = np.copy(y_temp), np.copy(t_temp)
    y_tran = frame_trans(y, t, fix_p = 1-mu_e,option='rota2init')
    y_m,t_m = primary_trans(y_tran, t, option='e2m')
    y_m_tran = frame_trans(y_m, t_m, fix_p = -mu_m, theta=theta0, \
                option='init2rota')
    return y_m_tran, t_m


def bcmfunc_em(y, t, theta_S0=0.):
    """
    地-月系下的考虑太阳摄动的方程
    """
    y1, y2, y3, y4, y5, y6 = y[0], y[1], y[2], y[3], y[4], y[5]

    epsilon = 1
    m_S = 328900.54 * epsilon
    a_S = 388.81114
    omega_S = 0.925195985520347
    mu_E = 1-mu_m
    mu_M = mu_m
    mu_S = m_S
    theta_S = -omega_S * t + theta_S0
    x_S = a_S * np.cos(theta_S)
    y_S = a_S * np.sin(theta_S)

    r_E = np.sqrt((y1+mu_M)**2+y2**2+y3**2)
    r_M = np.sqrt((y1-mu_E)**2+y2**2+y3**2)
    r_S = np.sqrt((y1-x_S)**2+(y2-y_S)**2+y3**2)
    c_E, c_M, c_S = mu_E/r_E**3, mu_M/r_M**3, mu_S/r_S**3
    alpha_S = m_S/a_S**3

    dy1 = y4
    dy2 = y5
    dy3 = y6
    dy4 = y1 + 2*y5 - c_E*(y1+mu_M) - c_M*(y1-mu_E) - c_S*(y1-x_S) - alpha_S*x_S
    dy5 = y2 - 2*y4 - c_E*y2 - c_M*y2 - c_S*(y2-y_S) - alpha_S*y_S
    dy6 = -c_E*y3 - c_M*y3 - c_S*y3

    return [dy1, dy2, dy3, dy4, dy5, dy6]

def bcmfunc_se(y, t, theta_M0=0.):
    """
    日-地/月系下的考虑月球摄动的方程
    """
    y1, y2, y3, y4, y5, y6 = y[0], y[1], y[2], y[3], y[4], y[5]

    epsilon = 1
    m_M = 3.733998734625702e-8 * epsilon
    a_M = 2.573565073532068e-3
    omega_M = 12.36886949284508
    mu_E = mu_e
    mu_M = m_M
    mu_S = 1-mu_e
    theta_M = omega_M * t + theta_M0
    x_M = a_M * np.cos(theta_M)
    y_M = a_M * np.sin(theta_M)

    r_E = np.sqrt((y1-mu_S)**2+y2**2+y3**2)
    r_M = np.sqrt((y1-x_M)**2+(y2-y_M)**2+y3**2)
    r_S = np.sqrt((y1+mu_E)**2+y2**2+y3**2)
    c_E, c_M, c_S = mu_E/r_E**3, mu_M/r_M**3, mu_S/r_S**3

    dy1 = y4
    dy2 = y5
    dy3 = y6
    dy4 = y1 + 2*y5 - c_S*(y1+mu_E) - c_E*(y1-mu_S) - c_M*(y1-x_M)
    dy5 = y2 - 2*y4 - c_S*y2 - c_E*y2 - c_M*(y2-y_M)
    dy6 = -c_S*y3 - c_E*y3 - c_M*y3

    return [dy1, dy2, dy3, dy4, dy5, dy6]

def bcmfunc_se_tran(t, y, theta_M0=0.):
    """
    日-地/月系下的考虑月球摄动的方程
    """
    y1, y2, y3, y4, y5, y6 = y[0], y[1], y[2], y[3], y[4], y[5]

    epsilon = 1
    m_M = 3.733998734625702e-8 * epsilon
    a_M = 2.573565073532068e-3
    omega_M = 12.36886949284508
    mu_E = mu_e
    mu_M = m_M
    mu_S = 1-mu_e
    theta_M = omega_M * t + theta_M0
    x_M = a_M * np.cos(theta_M)
    y_M = a_M * np.sin(theta_M)

    r_E = np.sqrt((y1-mu_S)**2+y2**2+y3**2)
    r_M = np.sqrt((y1-x_M)**2+(y2-y_M)**2+y3**2)
    r_S = np.sqrt((y1+mu_E)**2+y2**2+y3**2)
    c_E, c_M, c_S = mu_E/r_E**3, mu_M/r_M**3, mu_S/r_S**3

    dy1 = y4
    dy2 = y5
    dy3 = y6
    dy4 = y1 + 2*y5 - c_S*(y1+mu_E) - c_E*(y1-mu_S) - c_M*(y1-x_M)
    dy5 = y2 - 2*y4 - c_S*y2 - c_E*y2 - c_M*(y2-y_M)
    dy6 = -c_S*y3 - c_E*y3 - c_M*y3

    return [dy1, dy2, dy3, dy4, dy5, dy6]


def generate_amatrix_bcm_em(x, y, z, mu=mu_m, t=0., theta_S0=0.):
    """
    生成方程的系数矩阵A
    """
    epsilon = 1
    m_S = 328900.54 * epsilon
    a_S = 388.81114
    omega_S = 0.925195985520347
    mu_E = 1-mu_m
    mu_M = mu_m
    mu_S = m_S
    theta_S = -omega_S * t + theta_S0
    x_S = a_S * np.cos(theta_S)
    y_S = a_S * np.sin(theta_S)

    r_E = np.sqrt((x+mu_M)**2+y**2+z**2)
    r_M = np.sqrt((x-mu_E)**2+y**2+z**2)
    r_S = np.sqrt((x-x_S)**2+(y-y_S)**2+z**2)

    uxx = (1 - mu_E / r_E**3 - mu_M / r_M**3 + (3 * (x - mu_E)**2 *
           mu_M) / r_M**5 + (3 * (x + mu_M)**2 * mu_E) / r_E**5 -
           mu_S / r_S**3 + (3 * (x - x_S)**2 * mu_S) / r_S**5)

    uxy = ((3 * y * (x - mu_E) * mu_M) / r_M**5 + (3 * y * (x + mu_M) *
           mu_E) / r_E**5 + (3 * (y - y_S) * (x - x_S) * mu_S) / r_S**5)

    uxz = ((3 * z * (x - mu_E) * mu_M) / r_M**5 + (3 * z * (x + mu_M) *
           mu_E) / r_E**5 + (3 * z * (x - x_S) * mu_S) / r_S**5)

    uyy = (1 - mu_E / r_E**3 - mu_M / r_M**3 + (3 * y**2 *
           mu_M) / r_M**5 + (3 * y**2 * mu_E) / r_E**5 -
           mu_S / r_S**3 + (3 * (y - y_S)**2 * mu_S) / r_S**5)

    uyz = (3 * y * z * mu_E / r_E**5 +
           3 * y * z * mu_M / r_M**5 +
          (3 * (y - y_S) * z * mu_S) / r_S**5)

    uzz = (- mu_E / r_E**3 - mu_M / r_M**3 + (3 * z**2 * mu_E) / r_E**5 +
          (3 * z**2 * mu_M) / r_M**5 - mu_S / r_S**3 + (3 * z**2 * mu_S) /
          r_S**5)

    return np.array([[0, 0, 0, 1., 0, 0], [0, 0, 0, 0, 1., 0],
                     [0, 0, 0, 0, 0, 1.], [uxx, uxy, uxz, 0, 2., 0],
                     [uxy, uyy, uyz, -2., 0, 0], [uxz, uyz, uzz, 0, 0, 0]])


def dist2target(y, target='earth',model='se'):
    x, y, z = y[0], y[1], y[2]
    if target is 'earth':
        if model is 'se':
            pos = 1-mu_e
        elif model is 'em':
            pos = -mu_m
    elif target is 'moon':
        pos = 1-mu_m
    dist = np.sqrt((x - pos)**2 + y**2 + z**2)
    return dist

def xy2r(y, target='earth',model='se'):
    if target is 'earth':
        if model is 'se':
            center = np.array([[1-mu_e,0]])
        elif model is 'em':
            center = np.array([[-mu_m,0]])
    elif target is 'moon':
        center = np.array([[1-mu_m,0]])

    if len(y.shape) == 1:
        r_norm = np.linalg.norm(y[0:2]-center.reshape(2,))
        v_norm = np.linalg.norm(y[3:5])
        unit_vector = (y[0:2]-center.reshape(2,))/r_norm
        v_r = y[3:5].T.dot(unit_vector)
    else:
        r_norm = np.linalg.norm(y[0:2]-center.T,axis=0)
        v_norm = np.linalg.norm(y[3:5], axis=0)
        unit_vector = (y[0:2]-center.T)/r_norm
        v_r = y[3:5].T.dot(unit_vector).diagonal()
    v_t = np.sqrt(v_norm**2 - v_r**2)

    return r_norm, v_r, v_t

def plot_primaries(ax, option='em2d'):
    if option[0:2] == 'em':
        mu = mu_m
        L1 ,L2, L3 = L1_m, L2_m, L3_m
    elif option[0:2] == 'se':
        mu = mu_e
        L1 ,L2, L3 = L1_e, L2_e, L3_e
    if option[2:] == '2d':
        ax.plot([1-mu], [0], 'oy')
        ax.plot([-mu], [0], 'ob')
        ax.plot([L1,L2,L3], [0,0,0], 'xk')
    elif option[2:] == '3d':
        ax.plot([1-mu], [0], [0], 'oy')
        ax.plot([-mu], [0], [0], 'ob')
        ax.plot([L1,L2,L3], [0,0,0], [0,0,0], 'xk')

def two_body_energy(y, mu=mu_m):
    x, y, z, dotx, doty = y[0],y[1],y[2],y[3],y[4]
    pos = 1-mu_m
    v2 = np.sqrt((dotx-y)**2+(doty+x+mu-1)**2)
    r2 = np.sqrt((x - pos)**2 + y**2 + z**2)
    H2 = v2**2/2 - mu/r2
    #h2 = (x+mu-1)*(doty+x+mu-1)-y*(dotx-y)
    return H2, v2, r2

a = generate_amatrix(1, 1, 1, 0.001)
div = 1e-8
b = np.array([div,div,div,div,div,div])
c = np.dot(a,b)
print c
print np.linalg.norm(b),np.linalg.norm(c)