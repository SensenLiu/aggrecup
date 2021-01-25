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
import math

import platform


# --global variable--#
tf = 2
discretized_point_persecond = 50
pointnumber = tf * discretized_point_persecond  # 离散点数
controlstate_msg = controlstate()  # 要发布的控制量消息
# initialize
bvp_res = np.zeros((1, 7))
px_length = 0.0
pz_length = 0.0
vx_length = 0.0
vz_length = 0.0
va_length = 0.0
store_num = []
step = []
px_ini = 0.0
pz_ini = 0.0
vx_ini = 0.0
vz_ini = 0.0
va_ini = 0.0
pub = []
# --global variable--#

# TODO: interpolation? va_ini loop?
step = [0.1, 0.1, 0.1, 0.1, 0.5]  # discretized solving step ; type: list
interpolation = 1  # 0: round;   1: 2-ponits interpolation
model = 0  # 0: no va_ini loop ; 1: va_ini loop, air-resistance
boundary = np.zeros((2,5))


def pos_twist_callback(data):
    global px_ini, pz_ini, vx_ini, vz_ini, bvp_res, step, px_length, pz_length, vx_length, vz_length, va_length, store_num
    print("===callback===")
    px_ini = min(max(data.pose.pose.position.x, boundary[0, 0]), boundary[1, 0])  # relative pos, restrict
    pz_ini = min(max(data.pose.pose.position.z, boundary[0, 1]), boundary[1, 1])
    vx_ini = min(max(data.twist.twist.linear.x, boundary[0, 2]), boundary[1, 2])
    vz_ini = min(max(data.twist.twist.linear.z, boundary[0, 3]), boundary[1, 3])
    controlstate_msg.inicounter = 0
    controlstate_msg.discrepointpersecond = discretized_point_persecond
    controlstate_msg.arraylength = pointnumber
    controlstate_msg.thrustarray = []
    controlstate_msg.thetaarray = []

    start = time.time()

    # core code -- round the read list
    if interpolation == 0:
        px_index = round((px_ini - bvp_res[0, 0]) / step[0])
        pz_index = round((pz_ini - bvp_res[0, 1]) / step[1])
        vx_index = round((vx_ini - bvp_res[0, 2]) / step[2])
        vz_index = round((vz_ini - bvp_res[0, 3]) / step[3])
        va_index = round((va_ini - bvp_res[0, 4]) / step[4])
        if model == 0:  # no va_ini loop
            row = math.ceil(
                px_index * (pz_length * vx_length * vz_length) + pz_index * (vx_length * vz_length)
                + vx_index * vz_length + vz_index + 2.0)  # round up to an integer to avoid decreasing
        else:
            row = math.ceil(
                px_index * (pz_length * vx_length * vz_length * va_length) + pz_index * (
                        vx_length * vz_length * va_length)
                + vx_index * (
                            vz_length * va_length) + vz_index * va_length + va_index + 2.0)  # round up to an integer to avoid decreasing
        row = int(row)  # transform from float to int
        # print(row)
        controlstate_msg.thrustarray = bvp_res[row, 0:store_num]  # 0:5
        controlstate_msg.thetaarray = bvp_res[row, store_num:2 * store_num]  # 5:10
        controlstate_msg.stateXarray = bvp_res[row, 2 * store_num:3 * store_num]
        controlstate_msg.stateZarray = bvp_res[row, 3 * store_num:4 * store_num]
        controlstate_msg.stateVXarray = bvp_res[row, 4 * store_num:5 * store_num]
        controlstate_msg.stateVZarray = bvp_res[row, 5 * store_num:6 * store_num]
    else:
        px_index = [math.floor((px_ini - bvp_res[0, 0]) / step[0]), math.ceil((px_ini - bvp_res[0, 0]) / step[0])]
        pz_index = [math.floor((pz_ini - bvp_res[0, 1]) / step[1]), math.ceil((pz_ini - bvp_res[0, 1]) / step[1])]
        vx_index = [math.floor((vx_ini - bvp_res[0, 2]) / step[2]), math.ceil((vx_ini - bvp_res[0, 2]) / step[2])]
        vz_index = [math.floor((vz_ini - bvp_res[0, 3]) / step[3]), math.ceil((vz_ini - bvp_res[0, 3]) / step[3])]
        va_index = [math.floor((va_ini - bvp_res[0, 4]) / step[4]), math.ceil((va_ini - bvp_res[0, 4]) / step[4])]
        if model == 0:  # no va_ini loop
            row = [math.ceil(px_index[0] * (pz_length * vx_length * vz_length) + pz_index[0] * (vx_length * vz_length) +
                             vx_index[0] * vz_length + vz_index[0] + 2.0),
                   math.ceil(px_index[1] * (pz_length * vx_length * vz_length) + pz_index[1] * (vx_length * vz_length) +
                             vx_index[1] * vz_length + vz_index[1] + 2.0)]

        else:
            row = [math.ceil(px_index[0] * (pz_length * vx_length * vz_length * va_length) + pz_index[0] * (
                    vx_length * vz_length * va_length) + vx_index[0] * (vz_length * va_length) + vz_index[
                                 0] * va_length + va_index + 2.0),
                   math.ceil(px_index[1] * (pz_length * vx_length * vz_length * va_length) + pz_index[1] * (
                           vx_length * vz_length * va_length) + vx_index[1] * (vz_length * va_length) +
                             vz_index[1] * va_length + va_index + 2.0)]
        row = [int(row[0]), int(row[1])]  # transform from float to int
        temp_thrust = np.array([bvp_res[row[0], 0:store_num], bvp_res[row[1], 0:store_num]])  # (2,5)
        temp_theta = np.array([bvp_res[row[0], store_num:2 * store_num], bvp_res[row[1], store_num:2 * store_num]])
        temp_stateX = np.array(
            [bvp_res[row[0], 2 * store_num:3 * store_num], bvp_res[row[1], 2 * store_num:3 * store_num]])
        temp_stateZ = np.array(
            [bvp_res[row[0], 3 * store_num:4 * store_num], bvp_res[row[1], 3 * store_num:4 * store_num]])
        temp_stateVX = np.array(
            [bvp_res[row[0], 4 * store_num:5 * store_num], bvp_res[row[1], 4 * store_num:5 * store_num]])
        temp_stateVZ = np.array(
            [bvp_res[row[0], 5 * store_num:6 * store_num], bvp_res[row[1], 5 * store_num:6 * store_num]])

        # TODO: 求平均
        controlstate_msg.thrustarray = np.mean(temp_thrust, 0)  # mean by column
        controlstate_msg.thetaarray = np.mean(temp_theta, 0)
        controlstate_msg.stateXarray = np.mean(temp_stateX, 0)
        controlstate_msg.stateZarray = np.mean(temp_stateZ, 0)
        controlstate_msg.stateVXarray = np.mean(temp_stateVX, 0)
        controlstate_msg.stateVZarray = np.mean(temp_stateVZ, 0)
        print("thrust: %.5f" % controlstate_msg.thrustarray[0])
        print("pitch: %.5f" % controlstate_msg.thetaarray[0])
    # core code

    end = time.time()
    running_time = end - start
    print('time cost : %.5f sec' % running_time)

    # print(controlstate_msg.thetaarray[0])
    pub.publish(controlstate_msg)


def plane_vel_callback(data):
    global va_ini
    va_ini = min(max(data.twist.linear.x, boundary[0, 4]), boundary[1, 4])


def thread_offboard():
    global pub, discretized_point_persecond, bvp_res, step, px_length, pz_length, vx_length, vz_length, va_length, store_num, boundary
    rospy.init_node('bvp_control', anonymous=True)
    uav_id = rospy.get_param("~id", "")
    rate = rospy.Rate(100)

    rospy.Subscriber(uav_id + "current_relative_postwist",
                     Odometry, pos_twist_callback, queue_size=1)
    rospy.Subscriber(uav_id + "mavros/local_position/velocity_local",
                     TwistStamped, plane_vel_callback, queue_size=1)  # plane veocity

    pub = rospy.Publisher(
        uav_id + "bvp_controlstate", controlstate, queue_size=1)

    bvp_res = np.loadtxt(open("../data/multi_point/revise_bvp_table.csv", "rb"), delimiter=",", skiprows=0)
    print("===read===")
    # 读取列表
    store_num = int(bvp_res.shape[1] / 7)
    px_length = (bvp_res[1, 0] - bvp_res[0, 0]) / step[0]
    pz_length = (bvp_res[1, 1] - bvp_res[0, 1]) / step[1]
    vx_length = (bvp_res[1, 2] - bvp_res[0, 2]) / step[2]
    vz_length = (bvp_res[1, 3] - bvp_res[0, 3]) / step[3]
    va_length = (bvp_res[1, 4] - bvp_res[0, 4]) / step[4]
    boundary = np.array([bvp_res[0, 0:5], bvp_res[1, 0:5]])

    rospy.spin()


if __name__ == '__main__':  # 主函数shape
    thread_offboard()
