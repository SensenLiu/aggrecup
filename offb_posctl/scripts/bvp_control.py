#!/usr/bin/env python
# -*- coding: utf-8 -*-
# coding=utf-8
import sys

if sys.getdefaultencoding() != 'utf-8':
    reload(sys)
    sys.setdefaultencoding('utf-8')
import rospy

from nav_msgs.msg import Odometry
from geometry_msgs.msg import TwistStamped
from offb_posctl.msg import controlstate  # 发布自定义消息
from scipy.integrate import solve_bvp

import matplotlib.pyplot as plt
import numpy as np
import time

px_ini = -3.0
pz_ini = 0
vx_ini = 0.0
vz_ini = 0.0
va_ini = 0.0  # absolute velocity of plane
g = 9.8
k = np.array([50, 50])  # k[0]是x方向的系数，k[1]是z方向的系数
tf = 2  # 积分时间
discretized_point_persecond = 50
pointnumber = tf * discretized_point_persecond  # 离散点数
controlstate_msg = controlstate()  # 要发布的控制量消息
currentupdateflag = False  # 是否计算控制量
# c = np.array([1.5, 0.38])  # x和z方向的空气阻力系数
c = np.array([1.3, 0.25])  # x和z方向的空气阻力系数
model = 1  # 0: normal model ; 1:air-resistance model


def fun(t, y):
    dy = np.zeros(y.shape)

    u1 = -y[3] + g
    # u2 = (-g * y[2] - 2 * y[4] * y[5]) / (2 * y[4] ** 2 + 1)
    u2 = (-g * y[2] - 2 * k[1] * y[4] * y[5]) / (2 * k[0] * y[4] ** 2 + 1)

    # dy[0] = -2 * (y[4] * u2 + y[5]) * u2 - 2 * (y[4] + 3)
    dy[0] = -2 * k[1] * (y[4] * u2 + y[5]) * u2 - 2 * k[0] * (y[4] + 3)
    # dy[1] = -2 * (y[4] * u2 + y[5])
    dy[1] = -2 * k[1] * (y[4] * u2 + y[5])
    dy[2] = -y[0]
    dy[3] = -y[1]
    dy[4] = y[6]
    dy[5] = y[7]
    dy[6] = g * u2
    dy[7] = u1 - 9.8
    return dy


def bc(ya, yb):
    res = np.zeros(ya.shape)
    res[0] = ya[4] - px_ini
    res[1] = ya[5] - pz_ini
    res[2] = ya[6] - vx_ini
    res[3] = ya[7] - vz_ini
    res[4] = yb[0]
    res[5] = yb[1]
    res[6] = yb[2]
    res[7] = yb[3]
    return res


def do_process():
    t = np.linspace(0, tf, pointnumber)  # Define the initial mesh,行向量
    y_guess = np.zeros((8, t.size))  # initial guess for y
    res = solve_bvp(fun, bc, t, y_guess)
    t_plot = t
    y_plot = res.sol(t_plot)
    thrust = -y_plot[3] + g
    pitch = (-g * y_plot[2] - 2 * k[1] * y_plot[4] * y_plot[5]) / (2 * k[0] * y_plot[4] ** 2 + 1)
    return np.vstack((thrust, pitch, y_plot[4], y_plot[5], y_plot[6], y_plot[7]))


def pos_twist_callback(data):
    global px_ini, pz_ini, vx_ini, vz_ini, currentupdateflag
    px_ini = data.pose.pose.position.x  # relative pos
    pz_ini = data.pose.pose.position.z
    vx_ini = data.twist.twist.linear.x
    vz_ini = data.twist.twist.linear.z
    currentupdateflag = True


# ------------------------- take the air-resistance into consideration ----------------------------#
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
    return np.vstack((thrust, pitch, y_plot[5], y_plot[6], y_plot[7], y_plot[8]))


def plane_vel_callback(data):
    global va_ini
    va_ini = data.twist.linear.x


# ------------------------- take the air-resistance into consideration ----------------------------#


def thread_offboard():
    global currentupdateflag, discretized_point_persecond
    rospy.init_node('bvp_control', anonymous=True)
    uav_id = rospy.get_param("~id", "")
    rate = rospy.Rate(100)

    rospy.Subscriber(uav_id + "current_relative_postwist",
                     Odometry, pos_twist_callback)
    rospy.Subscriber(uav_id + "mavros/local_position/velocity_local",
                     TwistStamped, plane_vel_callback)  # plane veocity

    pub = rospy.Publisher(
        uav_id + "bvp_controlstate", controlstate, queue_size=10)

    while not (rospy.is_shutdown()):
        if currentupdateflag:
            controlstate_msg.inicounter = 0
            controlstate_msg.discrepointpersecond = discretized_point_persecond
            controlstate_msg.arraylength = pointnumber
            controlstate_msg.thrustarray = []
            controlstate_msg.thetaarray = []

            start = time.time()

            # core calculate code
            if model == 0:
                bvp_res = do_process()
            else:
                bvp_res = air_do_process()
            # core calculate code

            end = time.time()
            running_time = end - start
            print('time cost : %.5f sec' % running_time)

            # publish controlstate
            controlstate_msg.thrustarray = bvp_res[0, :]
            controlstate_msg.thetaarray = bvp_res[1, :]
            controlstate_msg.stateXarray = bvp_res[2, :]
            controlstate_msg.stateZarray = bvp_res[3, :]
            controlstate_msg.stateVXarray = bvp_res[4, :]
            controlstate_msg.stateVZarray = bvp_res[5, :]
            # print(controlstate_msg.thetaarray[0])
            pub.publish(controlstate_msg)
            currentupdateflag = False

        rate.sleep()


if __name__ == '__main__':  # 主函数
    thread_offboard()
