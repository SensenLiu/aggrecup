#!/usr/bin/env python
# -*- coding: utf-8 -*-
# coding=utf-8
import copy
import socket
import numpy as np
from scipy.optimize import minimize
import time
import datetime
import math
import matplotlib.pyplot as plt
import rospy

from nav_msgs.msg import Odometry
from geometry_msgs.msg import TwistStamped
from nav_msgs.msg import Path
from geometry_msgs.msg import PoseStamped
from offb_posctl.msg import controlstate  # 发布自定义消息

phi=1.57
normspeed=1.0
ay0=0
vy0=0
y0=0
az0=0
vz0=0
z0=0.5
aytf=-math.sin(phi)*9.8
vytf=normspeed*math.sin(phi)
ytf=2.0
aztf=math.cos(phi)*9.8-9.8
vztf=-normspeed*math.cos(phi)
ztf=2.0
meshpoint=np.linspace(1, 0.01, 5)
thrustmax=2*9.8
angleaccdmax=20
lbz=0.2
ubz=2.5
lbv=-5
ubv=5
currentupdateflag = False

# Objective
def J(x):
    return x[-1]

def fast_jac(x):
    jac = np.zeros_like(x)
    jac[-1]=1
    return jac

# Constraint
def eqmycon(x):
    global ay0, vy0, y0, az0, vz0, z0, aytf, vytf, ytf, aztf, vztf, ztf, meshpoint, thrustmax, angleaccdmax, lbz, lbv, ubv
    alpha_y=x[0]
    beta_y=x[1]
    gamma_y=x[2]

    alpha_z=x[3]
    beta_z=x[4]
    gamma_z=x[5]
    t=x[6]

    ceq1=alpha_y/6*t**3+beta_y/2*t**2+gamma_y*t+ay0-aytf
    ceq2=alpha_y/24*t**4+beta_y/6*t**3+gamma_y/2*t**2+ay0*t+vy0-vytf
    ceq3=alpha_y/120*t**5+beta_y/24*t**4+gamma_y/6*t**3+ay0/2*t**2+vy0*t+y0-ytf

    ceq4=alpha_z/6*t**3+beta_z/2*t**2+gamma_z*t+az0-aztf
    ceq5=alpha_z/24*t**4+beta_z/6*t**3+gamma_z/2*t**2+az0*t+vz0-vztf
    ceq6=alpha_z/120*t**5+beta_z/24*t**4+gamma_z/6*t**3+az0/2*t**2+vz0*t+z0-ztf
    return np.hstack((ceq1,ceq2,ceq3,ceq4,ceq5,ceq6)).ravel()

# Constraint
def ineqmycon(x):
    global ay0, vy0, y0, az0, vz0, z0, aytf, vytf, ytf, aztf, vztf, ztf, meshpoint, thrustmax, angleaccdmax, lbz, ubz, lbv, ubv
    alpha_y=x[0]
    beta_y=x[1]
    gamma_y=x[2]

    alpha_z=x[3]
    beta_z=x[4]
    gamma_z=x[5]
    t=x[6]
    tmesh=t*(np.array(meshpoint))
    angleacc=np.zeros_like(tmesh)
    for i in range(len(tmesh)):
        # print("i===",tmesh)
        t=tmesh[i]
        angleacc[i]=((((alpha_y*t**2)/2 + beta_y*t + gamma_y)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - (((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2)*((2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3))/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)**2 - ((beta_y + alpha_y*t)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - ((beta_z + alpha_z*t)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_z*t**2)/2 + beta_z*t + gamma_z))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)**2*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3)/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)

    t=x[6]
    thrust=np.sqrt(((alpha_y*tmesh**3)/6 + (beta_y*tmesh**2)/2 + gamma_y*tmesh + ay0)**2 + ((alpha_z*tmesh**3)/6 + (beta_z*tmesh**2)/2 + gamma_z*tmesh + az0 + 49/5)**2)
    c0=t
    # thrust constraints
    c1=2*9.8-thrust
    # print("c1----",c1.shape)
    # z's lower bound  constraints
    c2=-lbz+(alpha_z/120*tmesh**5+beta_z/24*tmesh**4+gamma_z/6*tmesh**3+az0/2*tmesh**2+vz0*tmesh+z0)
    c14=ubz-(alpha_z/120*tmesh**5+beta_z/24*tmesh**4+gamma_z/6*tmesh**3+az0/2*tmesh**2+vz0*tmesh+z0)

    # actuator constraints
    c3=angleacc*thrustmax/(4*angleaccdmax)-thrust/2+9.8
    c4=-angleacc*thrustmax/(4*angleaccdmax)+thrust/2
    c5=-angleacc*thrustmax/(4*angleaccdmax)-thrust/2+9.8
    c6=angleacc*thrustmax/(4*angleaccdmax)+thrust/2

    # phi belongs to [-1.57,1.57] constraints
    c7=((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 +9.8)
    c8=1
    c9=1
    if beta_z*beta_z>=2*alpha_z and alpha_z!=0 :
        t1=(-beta_z+math.sqrt(beta_z*beta_z-2*alpha_z))/alpha_z
        t2=(-beta_z-math.sqrt(beta_z*beta_z-2*alpha_z))/alpha_z
        if t1>=0 and t1<=t :
            c8=((alpha_z*t1**3)/6 + (beta_z*t1**2)/2 + gamma_z*t1 + az0 +9.8)
        if t2>=0 and t2<=t :
            c9=((alpha_z*t2**3)/6 + (beta_z*t2**2)/2 + gamma_z*t2 + az0 +9.8)
    #print('the value of t1 and t2 is',t1,t2)

    c10=-(alpha_y/24*tmesh**4+beta_y/6*tmesh**3+gamma_y/2*tmesh**2+ay0*tmesh+vy0-ubv)
    c11=-(lbv-(alpha_y/24*tmesh**4+beta_y/6*tmesh**3+gamma_y/2*tmesh**2+ay0*tmesh+vy0))
    c12=-(alpha_z/24*tmesh**4+beta_z/6*tmesh**3+gamma_z/2*tmesh**2+az0*tmesh+vz0-ubv)
    c13=-(lbv-(alpha_z/24*tmesh**4+beta_z/6*tmesh**3+gamma_z/2*tmesh**2+az0*tmesh+vz0))
    # print("--------", np.vstack((c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13)).shape)
    return np.hstack((c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14))

def pos_twist_callback(data):
    global vy0, y0, vz0, z0, currentupdateflag
    y0 = data.pose.pose.position.y  # relative pos
    z0 = data.pose.pose.position.z
    vy0 = data.twist.twist.linear.y
    vz0 = data.twist.twist.linear.z
    currentupdateflag = True

def plane_vel_callback(data):
    global va_ini
    va_ini = data.twist.linear.x

def droneImu_callback(data):
    global ay0, az0
    ay0 = data.twist.linear.y
    az0 = data.twist.linear.z

def main():
    global currentupdateflag
    constraint = [dict(type='eq', fun=eqmycon), dict(type='ineq', fun=ineqmycon)]
    Initial_guess=np.array([0,0,0,0,0,0,5])
    lb=-1000
    ub=1000
    mybounds=[(lb,ub),(lb,ub),(lb,ub),(lb,ub),(lb,ub),(lb,ub),(0,10)]

    controlfreq=30
    controlstate_msg = controlstate()  # 要发布的控制量消息
    planned_path=Path()
    planned_pose_stamped=PoseStamped()

    rospy.init_node('minimumsnap_control', anonymous=True)
    uav_id = rospy.get_param("~id", "")
    rate = rospy.Rate(100)

    rospy.Subscriber(uav_id + "current_relative_postwist",
                     Odometry, pos_twist_callback)
    # rospy.Subscriber(uav_id + "mavros/local_position/velocity_local",
    #                  TwistStamped, plane_vel_callback)  # plane veocity
    rospy.Subscriber(uav_id + "/mavros/imu/data",
                     TwistStamped, droneImu_callback)  # plane veocity

    pub = rospy.Publisher( uav_id + "ocp/control_state", controlstate, queue_size=1)
    path_pub = rospy.Publisher(uav_id + "plannedtrajectory",Path,queue_size=1)

    currentupdateflag=True
    while not (rospy.is_shutdown()):
        if currentupdateflag:
            start = time.time()
            result = minimize(J, Initial_guess, method='SLSQP', jac=fast_jac,tol=1e-4, bounds=mybounds,constraints=constraint)
            end = time.time()
            running_time = end - start
            if result.success:
                print('time cost : %.5f sec' % running_time)
                Initial_guess=result.x
                controlstate_msg.inicounter = 10
                controlstate_msg.discrepointpersecond = controlfreq
                controlstate_msg.arraylength = round(result.x[-1]*controlfreq)

                times=np.linspace(0, 1, controlstate_msg.arraylength)*result.x[-1]
                alpha_y=result.x[0]
                beta_y=result.x[1]
                gamma_y=result.x[2]

                alpha_z=result.x[3]
                beta_z=result.x[4]
                gamma_z=result.x[5]

                y=alpha_y/120*times**5+beta_y/24*times**4+gamma_y/6*times**3+ay0/2*times**2+vy0*times+y0
                vy=alpha_y/24*times**4+beta_y/6*times**3+gamma_y/2*times**2+ay0*times+vy0
                ay=alpha_y/6*times**3+beta_y/2*times**2+gamma_y*times+ay0

                z=alpha_z/120*times**5+beta_z/24*times**4+gamma_z/6*times**3+az0/2*times**2+vz0*times+z0
                vz=alpha_z/24*times**4+beta_z/6*times**3+gamma_z/2*times**2+az0*times+vz0
                az=alpha_z/6*times**3+beta_z/2*times**2+gamma_z*times+az0

                controlstate_msg.stateXarray = np.zeros_like(times)
                controlstate_msg.stateYarray = y
                controlstate_msg.stateZarray = z
                controlstate_msg.stateVXarray = np.zeros_like(times)
                controlstate_msg.stateVYarray = vy
                controlstate_msg.stateVZarray = vz
                controlstate_msg.stateAXarray = np.zeros_like(times)
                controlstate_msg.stateAYarray = ay
                controlstate_msg.stateAZarray = az

                pub.publish(controlstate_msg)
                currentupdateflag = False

                planned_path.header.stamp=rospy.Time.now()
                planned_path.header.frame_id="ground_link"
                planned_path.poses=[]
                for i in range(0, int(controlstate_msg.arraylength)):
                    planned_pose_stamped.pose.position.x=controlstate_msg.stateXarray[i]
                    planned_pose_stamped.pose.position.y=controlstate_msg.stateYarray[i]
                    planned_pose_stamped.pose.position.z=controlstate_msg.stateZarray[i]
                    planned_pose_stamped.header.stamp=rospy.Time.now()
                    # planned_pose_stamped.header.frame_id="ground_link"
                    planned_pose_stamped.header.seq=i
                    planned_path.poses.append(copy.deepcopy(planned_pose_stamped))
                path_pub.publish(planned_path)


    rate.sleep()
    # times=np.linspace(0,1,100)*result.x[-1]
    #
    # alpha_y=result.x[0]
    # beta_y=result.x[1]
    # gamma_y=result.x[2]
    #
    # alpha_z=result.x[3]
    # beta_z=result.x[4]
    # gamma_z=result.x[5]
    #
    # y=alpha_y/120*times**5+beta_y/24*times**4+gamma_y/6*times**3+ay0/2*times**2+vy0*times+y0
    # vy=alpha_y/24*times**4+beta_y/6*times**3+gamma_y/2*times**2+ay0*times+vy0
    # ay=alpha_y/6*times**3+beta_y/2*times**2+gamma_y*times+ay0
    #
    #
    # z=alpha_z/120*times**5+beta_z/24*times**4+gamma_z/6*times**3+az0/2*times**2+vz0*times+z0
    # vz=alpha_z/24*times**4+beta_z/6*times**3+gamma_z/2*times**2+az0*times+vz0
    # az=alpha_z/6*times**3+beta_z/2*times**2+gamma_z*times+az0
    # a=np.sqrt(az**2+ay**2)
    # thurst=np.sqrt((az+9.8)**2+ay**2)
    # phiseries=-np.arctan(ay/(az+9.8))
    # # print("az--------", az)
    # angleacc=np.zeros_like(times)
    # for i in range(len(times)):
    #     t=times[i]
    #     angleacc[i]=((((alpha_y*t**2)/2 + beta_y*t + gamma_y)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - (((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2)*((2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3))/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)**2 - ((beta_y + alpha_y*t)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - ((beta_z + alpha_z*t)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_z*t**2)/2 + beta_z*t + gamma_z))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)**2*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3)/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)
    #
    # F1=angleacc*thrustmax/(4*angleaccdmax)+thurst/2
    # F2=-angleacc*thrustmax/(4*angleaccdmax)+thurst/2
    #
    # plotlinewidth = 2
    # plotfontsize = 16
    # plt.subplot(2,2,1)
    # plt.plot(times,y, color='blue',LineWidth=plotlinewidth,label="y")
    # plt.plot(times,vy, color='green',LineWidth=plotlinewidth,label="vy")
    # plt.plot(times,ay, color='black', LineWidth=plotlinewidth,label="ay")
    # plt.plot(times,phiseries, color='yellow',LineWidth=plotlinewidth,label="phi")
    # plt.legend(loc="best")
    #
    # plt.subplot(2,2,2)
    # plt.plot(times,z, color='blue',LineWidth=plotlinewidth,label="z")
    # plt.plot(times,vz, color='green',LineWidth=plotlinewidth,label="vz")
    # plt.plot(times,az, color='black', LineWidth=plotlinewidth,label="az")
    # plt.plot(times,thurst, color='yellow',LineWidth=plotlinewidth,label="thurst")
    # plt.legend(loc="best")
    #
    # plt.subplot(2,2,3)
    # plt.plot(-y,z, color='blue',LineWidth=plotlinewidth,label="y-z")
    # plt.legend(loc="best")
    #
    # plt.subplot(2,2,4)
    # plt.plot(times,F1, color='blue',LineWidth=plotlinewidth,label="F1")
    # plt.plot(times,F2, color='black', LineWidth=plotlinewidth,label="F2")
    # plt.legend(loc="best")
    #
    # # print(res)
    # # core calculate code
    # plt.show()

if __name__ == '__main__':  # 主函数
    main()
