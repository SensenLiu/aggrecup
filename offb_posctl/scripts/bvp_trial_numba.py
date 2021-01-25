#!/usr/bin/env python
# -*- coding: utf-8 -*-
# coding=utf-8

from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
import numpy as np
import time
from numba import jit, float64
class bvp_class():
    px_ini = -3.23815
    pz_ini = 0.3406
    vx_ini = -0.17139
    vz_ini = 0.18496
    va_ini = -0.16429   # absolute velocity of plane
    g = 9.8
    k = np.array([50, 50])
    tf = 3  # 积分时间
    pointnumber = tf * 50  # 离散点数
    # c = np.array([1.5, 0.38])  # x和z方向的空气阻力系数
    c = np.array([0, 0])  # x和z方向的空气阻力系数
    model = 1  # 0: normal model ; 1:air-resistance model

    def _init(self):
        self.px_ini=-3

    # def fun(t, y):
    #     dy = np.zeros(y.shape)
    #
    #     u1 = -y[3] + g
    #     # u2 = (-g * y[2] - 2 * y[4] * y[5]) / (2 * y[4] ** 2 + 1)
    #     u2 = (-g * y[2] - 2 * k[1] * y[4] * y[5]) / (2 * k[0] * y[4] ** 2 + 1)
    #
    #     # dy[0] = -2 * (y[4] * u2 + y[5]) * u2 - 2 * (y[4] + 3)
    #     dy[0] = -2 * k[1] * (y[4] * u2 + y[5]) * u2 - 2 * k[0] * (y[4] + 3)
    #     # dy[1] = -2 * (y[4] * u2 + y[5])
    #     dy[1] = -2 * k[1] * (y[4] * u2 + y[5])
    #     dy[2] = -y[0]
    #     dy[3] = -y[1]
    #     dy[4] = y[6]
    #     dy[5] = y[7]
    #     dy[6] = g * u2
    #     dy[7] = u1 - 9.8
    #     return dy
    #
    #
    # def bc(ya, yb):
    #     res = np.zeros(ya.shape)
    #
    #     res[0] = ya[4] - px_ini
    #     res[1] = ya[5] - pz_ini
    #     res[2] = ya[6] - vx_ini
    #     res[3] = ya[7] - vz_ini
    #     res[4] = yb[0]
    #     res[5] = yb[1]
    #     res[6] = yb[2]
    #     res[7] = yb[3]
    #     return res
    #
    #
    # def do_process():
    #     t = np.linspace(0, tf, pointnumber)
    #     y_guess = np.zeros((8, t.size))  # initial guess for y
    #     res = solve_bvp(fun, bc, t, y_guess)
    #     t_plot = t
    #     y_plot = res.sol(t_plot)
    #     thrust = -y_plot[3] + g
    #     # pitch = (-g * y_plot[2] - 2 * y_plot[4] * y_plot[5]) / (2 * y_plot[4] ** 2 + 1)
    #     pitch = (-g * y_plot[2] - 2 * k[1] * y_plot[4] * y_plot[5]) / (2 * k[0] * y_plot[4] ** 2 + 1)
    #     return np.vstack((thrust, pitch, y_plot[4, :], y_plot[5, :], y_plot[6, :], y_plot[7, :]))


    # ------------------------- take the air-resistance into consideration ----------------------------#
    @jit(float64[:, :](float64[:], float64[:, :]), nopython=True)
    def air_fun(t, y):
        dy = np.zeros(y.shape)
        u1 = -y[3] + g
        u2 = (-g * (y[2] + y[4]) - 2 * k[1] * y[5] * y[6]) / (2 * k[0] * y[5] ** 2 + 1)
        dy[0] = -2 * k[1] * (y[5] * u2 + y[6]) * u2 - 2 * k[0] * (y[5] + 3)  # lamda1_dot
        dy[1] = -2 * k[1] * (y[5] * u2 + y[6])  # lamda2_dot
        dy[2] = -y[0]  # lamda3_dot
        dy[3] = -y[1]  # lamda4_dot
        dy[4] = c[0] * y[2] + c[1] * y[3] + c[0] * y[4]  # lamda5_dot
        dy[5] = y[7]  # x1_dot
        dy[6] = y[8]  # x2_dot
        dy[7] = g * u2 - c[0] * y[9]  # x3_dot
        dy[8] = u1 - 9.8 - c[1] * y[9]  # x4_dot
        dy[9] = g * u2 - c[0] * y[9]  # x5_dot
        return dy


    @jit(float64[:](float64[:], float64[:]), nopython=True)
    def air_bc(ya, yb):
        res = np.zeros(ya.shape)

        res[0] = ya[5] - px_ini
        res[1] = ya[6] - pz_ini
        res[2] = ya[7] - vx_ini
        res[3] = ya[8] - vz_ini
        res[4] = ya[9] - va_ini
        res[5] = yb[0]
        res[6] = yb[1]
        res[7] = yb[2]
        res[8] = yb[3]
        res[9] = yb[4]
        return res


    def air_do_process():
        t = np.linspace(0, tf, pointnumber)  # Define the initial mesh,行向量
        y_guess = np.zeros((10, t.size))  # initial guess for y
        res = solve_bvp(air_fun, air_bc, t, y_guess)
        t_plot = t
        y_plot = res.sol(t_plot)
        thrust = -y_plot[3] + g
        pitch = (-g * y_plot[2] - 2 * k[1] * y_plot[5] * y_plot[6]) / (2 * k[0] * y_plot[5] ** 2 + 1)
        return np.vstack((thrust, pitch, y_plot[5], y_plot[6], y_plot[7], y_plot[8], y_plot[9]))


# ------------------------- take the air-resistance into consideration ----------------------------#
a=bvp_class(3)
t_plot = np.linspace(0, a.tf, a.pointnumber)
start = time.time()
# core calculate code
for i in range(100):
    if a.model == 0:
        bvp_res = a.do_process()
    else:
        bvp_res = a.air_do_process()
# core calculate code
end = time.time()
running_time = end - start
print('time cost : %.5f sec' % running_time)
# print(bvp_res[6])

plt.figure(1)
plt.subplot(2, 2, 1)
plt.plot(t_plot, bvp_res[2], 'r--', t_plot, bvp_res[3], 'b-', linewidth=2)
plt.legend(labels=['xp-xt', 'zp-zt'])
plt.xlabel("t")
plt.grid()
plt.subplot(2, 2, 2)
plt.plot(t_plot, bvp_res[4], 'r--', t_plot, bvp_res[5], 'b-', linewidth=2)
plt.legend(labels=['vpx-vtx', 'vpz - vcz'])
plt.xlabel("t")
plt.grid()
plt.subplot(2, 2, 3)
plt.plot(t_plot, bvp_res[0], 'r--', linewidth=2)
plt.legend(labels=['thrust'])
plt.xlabel("t")
plt.grid()
plt.subplot(2, 2, 4)
plt.plot(t_plot, bvp_res[1], 'r--', linewidth=2)
plt.legend(labels=['pitch'])
plt.xlabel("t")
plt.grid()
plt.show()
# plt.figure(2)
# plt.plot(t_plot, bvp_res[2], 'r--', t_plot, bvp_res[3], 'b-', linewidth=2)
# plt.legend(labels=['px-cx', 'pz-z0'])
# plt.xlabel("t")
# plt.grid()
# plt.show()