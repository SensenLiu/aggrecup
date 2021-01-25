#!/usr/bin/env python
# -*- coding: utf-8 -*-
# coding=utf-8
import socket
from scipy.optimize import minimize
import numpy as np
import time
from numba import jit, float64
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import datetime

# global variable start >>>
n = 7
t0 = 0
tf = 2
discretized_point_persecond = 50
pointnumber = tf * discretized_point_persecond  # 离散点数
currentupdateflag = False  # 是否计算控制量
k = np.array([50, 50])
c = np.array([0, 0])  # air drag effect in x & z
co = 0.5 * (tf - t0)
g = 9.8
px_ini = -3
pz_ini = 0
vx_ini = 0
vz_ini = 0
va_ini = 0  # absolute velocity of plane
# ini = np.array([[px_ini], [pz_ini], [vx_ini], [vz_ini], [va_ini]])
state_get_flag = False
# global variable start >>>

# D matrix
D = np.loadtxt(open("../data/D.csv", "rb"), delimiter=",", skiprows=0)  # array
#  Gauss weights
omega = np.loadtxt(open("../data/omega.csv", "rb"), delimiter=",", skiprows=0)  # array
# Lagrange coefficient of x
L1 = np.loadtxt(open("../data/L1.csv", "rb"), delimiter=",", skiprows=0)  # array
# Lagrange coefficient of u
L2 = np.loadtxt(open("../data/L2.csv", "rb"), delimiter=",", skiprows=0)  # array



# Objective
@jit(float64(float64[:]), nopython=True)
def J(x):
    X1 = x[0: n]
    X2 = x[n: 2 * n]
    U1 = x[5 * n: 6 * n]
    U2 = x[6 * n: 7 * n]
    return co * 0.5 * np.dot(omega, (
            0.5 * (U1 - 9.8) ** 2 + 0.5 * U2 ** 2 + k[0] * (X1 + 3) ** 2 + k[1] * (X1 * U2 + X2) ** 2))

# the derivative of objective function J
@jit(float64[:](float64[:]), nopython=True)
def fast_jac(x):
    h = 1e-11
    N = x.shape[0]
    jac = np.zeros_like(x)
    f_0 = J(x)
    for i in range(N):
        x_d = np.copy(x)
        x_d[i] += h
        f_d = J(x_d)
        jac[i] = (f_d - f_0) / h
    return jac


# Constraint
@jit(float64[:](float64[:]), nopython=True)
def mycon(x):
    global px_ini, pz_ini, vx_ini, vz_ini, va_ini
    X1 = x[0: n]
    X2 = x[n: 2 * n]
    X3 = x[2 * n: 3 * n]
    X4 = x[3 * n: 4 * n]
    X5 = x[4 * n: 5 * n]
    U1 = x[5 * n: 6 * n]
    U2 = x[6 * n: 7 * n]
    print('===================????', px_ini)
    Ceq1 = np.dot(D, np.append(px_ini, X1)) - co * X3
    Ceq2 = np.dot(D, np.append(pz_ini, X2)) - co * X4
    Ceq3 = np.dot(D, np.append(vx_ini, X3)) - co * (g * U2 - c[0] * X5)
    Ceq4 = np.dot(D, np.append(vz_ini, X4)) - co * (U1 - g - c[1] * X5)
    Ceq5 = np.dot(D, np.append(va_ini, X5)) - co * (g * U2 - c[0] * X5)
    return np.hstack((Ceq1, Ceq2, Ceq3, Ceq4, Ceq5))


def do_process(result):
    global tau
    x = result.x.reshape(7, n)
    print('===================!!!!', px_ini)
    ini = np.array([[px_ini], [pz_ini], [vx_ini], [vz_ini], [va_ini]])
    # print('ini.{}'.format(ini))
    poly_x = np.dot(np.hstack((ini, x[0:5, :])), L1)  # 拟合出的x的系数矩阵
    poly_u = np.dot(x[5:7, :], L2)  # 拟合出的u的系数矩阵

    # 将数据代入系数矩阵求x和u
    x1 = np.polyval(poly_x[0], tau)
    x2 = np.polyval(poly_x[1], tau)
    x3 = np.polyval(poly_x[2], tau)
    x4 = np.polyval(poly_x[3], tau)
    x5 = np.polyval(poly_x[4], tau)
    u1 = np.polyval(poly_u[0], tau)
    u2 = np.polyval(poly_u[1], tau)

    return np.vstack((x1, x2, x3, x4, u1, u2))


def parse(data):  # 解析收到的client的px等数据
    global state_get_flag
    if len(data) > 6:  # 判断是否包含至少包头
        for i in range(len(data)):
            if data[i:i + 3].decode() == 'LEN':
                Length = int(data[i + 3:i + 6].decode())
                # print('data:{}'.format(data))
                # print('time now:{}'.format(time.time()))
                if len(data[i:]) >= (Length + 6):  # 消息包含包头+state
                    msg = eval(data[i + 6:i + 6 + Length].decode())  # 直到处理完,最新msg
                    print('msg:{}'.format(msg))
                    if len(msg) == 6:
                        state_get_flag = True
                    if len(data[i + 6+Length:]) < Length + 6:  # 剩下的不够一条消息的长度
                        break
                else:
                    break
        try:
            return data[Length + i + 6:], msg  # 返回剩余不能构成一帧数据的data
        except:
            print('----data:{}----'.format(data))
            return b''
    else:
        return b''

    pass


def main():
    global currentupdateflag, discretized_point_persecond, tau, state_get_flag, msg
    global px_ini, pz_ini, vx_ini, vz_ini, va_ini

    constraint = [dict(type='eq', fun=mycon)]
    tau = np.linspace(-1, 1, pointnumber)
    # t = 0.5 * (tf - t0) * tau + 0.5 * (tf + t0)
    # while True:
    state_get_flag = True
    if state_get_flag:

        px_ini = -3.5
        pz_ini = 0.0
        vx_ini = -0.1
        vz_ini = 0
        va_ini = 0
        print('px_ini:{}; pz_ini:{}; vx_ini:{}; vz_ini:{}; va_ini:{};'.format(px_ini, pz_ini, vx_ini, vz_ini, va_ini))

        start = time.time()
        # core calculate code
        result = minimize(J, np.zeros((7 * n)), method='SLSQP', tol=1e-4, constraints=constraint, jac=fast_jac)
        print(result)
        res = do_process(result)
        # print(res)
        # core calculate code
        end = time.time()
        running_time = end - start
        print('time cost : %.5f sec' % running_time)

        ## core part 1 >>>>
        time_now = time.time()
        thrust_pitch_x1234 = [res[4, 0:20].tolist(), res[5, 0:20].tolist(), res[0, 0:20].tolist(), res[1, 0:20].tolist(),
                              res[2, 0:20].tolist(), res[3, 0:20].tolist(), time_now]
        plt.plot(tau*(tf-t0)/2.0+(tf+t0), res[0, 0:100])
        plt.show()


if __name__ == '__main__':  # 主函数
    main()
