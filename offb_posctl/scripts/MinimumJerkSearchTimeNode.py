#!/usr/bin/env python
# -*- coding: utf-8 -*-
# coding=utf-8
from __future__ import division
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
from sensor_msgs.msg import Imu
from nav_msgs.msg import Path
from geometry_msgs.msg import PoseStamped
from offb_posctl.msg import controlstate
from offb_posctl.msg import wallstate
from pyquaternion import Quaternion
# print(offb_posctl.__file__)

phi=1.57
normspeed=-0.15
tangentialspeed=0.8
ytf_wall=2.00
ztf_wall=2.00 # true is 1.94
vytf_wall=0.0
vztf_wall=0.0
cuplength=0.1

# parabolictime=phi/25.0
parabolictime=0
print(parabolictime*30)
approachthrust=9.8
a_normal=approachthrust-9.8*math.cos(phi)
a_tangential=-9.8*math.sin(phi)
normspeed=normspeed-a_normal*parabolictime
tangentialspeed=tangentialspeed-a_tangential*parabolictime
normaloff=-(normspeed*parabolictime+0.5*a_normal*parabolictime**2)+cuplength
tangentialoff=-(tangentialspeed*parabolictime+0.5*a_tangential*parabolictime**2)


ay0=0
vy0=0
y0=0
az0=0
vz0=0
z0=0.5
# ay0=-1.979504701384323
# vy0=3.7674505710601807
# y0=0.16079023480415344
# az0=-6.541593055514754
# vz0=3.7674505710601807
# z0=1.8277690410614014

delta_vytf=-normspeed*math.sin(phi)+tangentialspeed*math.cos(phi)
delta_ytf=-normaloff*math.sin(phi)+tangentialoff*math.cos(phi)
delta_vztf=normspeed*math.cos(phi)+tangentialspeed*math.sin(phi)
delta_ztf=normaloff*math.cos(phi)+tangentialoff*math.sin(phi)

aytf=-math.sin(phi)*9.8
vytf=delta_vytf+0
ytf=ytf_wall+delta_ytf

aztf=math.cos(phi)*9.8-9.8
vztf=delta_vztf+0
ztf=ztf_wall+delta_ztf
meshpoint=np.linspace(1, 0.01, 20)
thrustmax=2*9.8
angleaccdmax=25
lbz=0.2
ubz=2.5
lbv=-5
ubv=5
currentupdateflag = False
wallstateupdateflag = False
Init_guess=0.1
increaseratio=1.5

wallstate_msg = wallstate()  # 要发布的控制量消息
lastwallupdatetime=0.0
lastdroneupdatetime=0.0

lastsolved_time=0.0
lastsolveduration=0.0
# Constraint
def ineqmycon(x):
    global ay0, vy0, y0, az0, vz0, z0, aytf, vytf, ytf, aztf, vztf, ztf, delta_vytf, delta_ytf, delta_vztf, delta_ztf, ytf_wall,ztf_wall,vytf_wall,vztf_wall,meshpoint, thrustmax, angleaccdmax, lbz, ubz, lbv, ubv,lastwallupdatetime, wallstate_msg
    t=x
    # if((time.time()-lastwallupdatetime)*wallstate_msg.discrepointpersecond>=wallstate_msg.arraylength):
    #     lastwallupdatetime=False
    #     return False
    # wallindex=int(min(math.floor((t+time.time()-lastwallupdatetime)*wallstate_msg.discrepointpersecond),wallstate_msg.arraylength-1))
    # # print(wallindex)
    # if wallindex>=wallstate_msg.arraylength:
    #     return False
    # ytf_wall=wallstate_msg.stateYarray[wallindex]
    # vytf_wall=wallstate_msg.stateVYarray[wallindex]
    # ztf_wall=wallstate_msg.stateZarray[wallindex]
    # vztf_wall=wallstate_msg.stateVZarray[wallindex]

    # ytf=ytf_wall+delta_ytf
    # vytf=vytf_wall+delta_vytf
    # ztf=ztf_wall+delta_ztf
    # vztf=vztf_wall+delta_vztf

    # print("t----wallindex----ytf:",t,wallindex,ytf)
    tarray=np.array([[60/t**3,-360/t**4,720/t**5],[-24/t**2,168/t**3,-360/t**4],[3/t,-24/t**2,60/t**3]])

    alpha_y,beta_y,gamma_y=np.dot(tarray,np.array([aytf-ay0,vytf-vy0-ay0*t,ytf-y0-vy0*t-0.5*ay0*t**2]))
    alpha_z,beta_z,gamma_z=np.dot(tarray,np.array([aztf-az0,vztf-vz0-az0*t,ztf-z0-vz0*t-0.5*az0*t**2]))

    tmesh=t*(np.array(meshpoint))
    angleacc=np.zeros_like(tmesh)
    for i in range(len(tmesh)):
        # print("i===",tmesh)
        t=tmesh[i]
        angleacc[i]=((((alpha_y*t**2)/2 + beta_y*t + gamma_y)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - (((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2)*((2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3))/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)**2 - ((beta_y + alpha_y*t)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - ((beta_z + alpha_z*t)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_z*t**2)/2 + beta_z*t + gamma_z))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)**2*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3)/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)

    thrust=np.sqrt(((alpha_y*tmesh**3)/6 + (beta_y*tmesh**2)/2 + gamma_y*tmesh + ay0)**2 + ((alpha_z*tmesh**3)/6 + (beta_z*tmesh**2)/2 + gamma_z*tmesh + az0 + 49/5)**2)
    # thrust constraints
    # c1=2*9.8-thrust
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
    # c7=((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 +9.8)
    # c8=1
    # c9=1
    # if beta_z*beta_z>=2*alpha_z and alpha_z!=0 :
    #     t1=(-beta_z+math.sqrt(beta_z*beta_z-2*alpha_z))/alpha_z
    #     t2=(-beta_z-math.sqrt(beta_z*beta_z-2*alpha_z))/alpha_z
    #     if t1>=0 and t1<=t :
    #         c8=((alpha_z*t1**3)/6 + (beta_z*t1**2)/2 + gamma_z*t1 + az0 +9.8)
    #     if t2>=0 and t2<=t :
    #         c9=((alpha_z*t2**3)/6 + (beta_z*t2**2)/2 + gamma_z*t2 + az0 +9.8)
    #print('the value of t1 and t2 is',t1,t2)

    c10=-(alpha_y/24*tmesh**4+beta_y/6*tmesh**3+gamma_y/2*tmesh**2+ay0*tmesh+vy0-ubv)
    c11=-(lbv-(alpha_y/24*tmesh**4+beta_y/6*tmesh**3+gamma_y/2*tmesh**2+ay0*tmesh+vy0))
    c12=-(alpha_z/24*tmesh**4+beta_z/6*tmesh**3+gamma_z/2*tmesh**2+az0*tmesh+vz0-ubv)
    c13=-(lbv-(alpha_z/24*tmesh**4+beta_z/6*tmesh**3+gamma_z/2*tmesh**2+az0*tmesh+vz0))
    # print("--------t,flag", t,(np.hstack((c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14))))
    return (np.hstack((c2,c3,c4,c5,c6,c10,c11,c12,c13,c14))>-0.05).all()

def decreasedtimemethod():
    global lastsolved_time,lastsolveduration,lastwallupdatetime,wallstate_msg,vytf, ytf, vztf, ztf, delta_vytf, delta_ytf, delta_vztf, delta_ztf, ytf_wall,ztf_wall,vytf_wall,vztf_wall
    # rendezvousindex=int(min(math.floor((lastsolved_time+lastsolveduration-lastwallupdatetime)*wallstate_msg.discrepointpersecond),wallstate_msg.arraylength-1))
    #
    # ytf_wall=wallstate_msg.stateYarray[rendezvousindex]
    # vytf_wall=wallstate_msg.stateVYarray[rendezvousindex]
    # ztf_wall=wallstate_msg.stateZarray[rendezvousindex]
    # vztf_wall=wallstate_msg.stateVZarray[rendezvousindex]
    #
    # ytf=ytf_wall+delta_ytf
    # vytf=vytf_wall+delta_vytf
    # ztf=ztf_wall+delta_ztf
    # vztf=vztf_wall+delta_vztf

def pos_twist_callback(data):
    global vy0, y0, vz0, z0, ay0,az0,currentupdateflag,targetstartmoveflag
    y0 = data.pose.pose.position.y  # relative pos
    z0 = data.pose.pose.position.z
    vy0 = data.twist.twist.linear.y
    vz0 = data.twist.twist.linear.z
    ay0=data.twist.twist.angular.y
    az0=data.twist.twist.angular.z
    # print(ay0)
    currentupdateflag = True


def droneImu_callback(data):
    global ay0, az0,lastdroneupdatetime
    # ay0 = data.linear_acceleration.y
    # az0 = data.linear_acceleration.z-9.8
    qcurrent = Quaternion(data.orientation.w, data.orientation.x, data.orientation.y, data.orientation.z)
    qaccinworld=(qcurrent).rotate(Quaternion(vector=[data.linear_acceleration.x,data.linear_acceleration.y,data.linear_acceleration.z]))
    # ay0=0.2*min(max(qaccinworld.y,-2*9.8),2*9.8)+0.8*ay0
    # az0=0.2*min(max(qaccinworld.z-9.8,-9.8),9.8)+0.8*az0
    # ay0=0
    # az0=0
    # print ("ay0",ay0)
    lastdroneupdatetime=time.time()

def wallstate_callback(data):
    global wallstate_msg,wallstateupdateflag,lastwallupdatetime
    wallstate_msg=data
    wallstateupdateflag=True
    lastwallupdatetime=time.time()
    # print ("ay0",ay0)

def main():
    global currentupdateflag, wallstateupdateflag, Init_guess,increaseratio, targetstartmoveflag, ytf, ztf,lastsolved_time,lastdroneupdatetime,lastsolveduration
    controlfreq=30
    controlstate_msg = controlstate()  # 要发布的控制量消息
    planned_path=Path()
    target_path=Path()
    planned_pose_stamped=PoseStamped()
    target_pose_stamped=PoseStamped()
    recv_acc=Imu()

    rospy.init_node('minimumsnap_control', anonymous=True)
    uav_id = rospy.get_param("~id", "")
    rate = rospy.Rate(100)

    rospy.Subscriber(uav_id + "current_relative_postwist",
                     Odometry, pos_twist_callback)
    # rospy.Subscriber(uav_id + "mavros/local_position/velocity_local",
    #                  TwistStamped, plane_vel_callback)  # plane veocity
    rospy.Subscriber(uav_id + "/mavros/imu/data",
                     Imu, droneImu_callback)  # plane veocity
    rospy.Subscriber(uav_id + "predictedwall_array",
                     wallstate, wallstate_callback)  # plane veocity

    pub = rospy.Publisher( uav_id + "ocp/control_state", controlstate, queue_size=1)
    path_pub = rospy.Publisher(uav_id + "plannedtrajectory",Path,queue_size=1)
    acc_pub = rospy.Publisher(uav_id +"/acc_data_inworld",Imu,queue_size=1)

    currentupdateflag=True
    searchstep,leftnode,rightnode,solveflag=1, 0.1, 10, False

    while not (rospy.is_shutdown()):
        # if targetstartmoveflag:
        #     ytf=ytf+(targetvelamplitude*math.sin(2*math.pi*2*targetmovecounter*100)+targetbasevelocity)*0.01
        #     targetmovecounter=targetmovecounter+1
        # if currentupdateflag and wallstateupdateflag:
        if currentupdateflag :
            Init_guess=leftnode*0.5
            leftnode=Init_guess
            rightnode=min(rightnode*increaseratio,10)
            searchstep=max((rightnode-leftnode)/5, 0.01)
            solveflag=False

            while rightnode-leftnode>0.2 or solveflag==False:
                solveflag=ineqmycon(leftnode)
                # print("solveflag,leftnode",solveflag,leftnode)
                if solveflag==True:
                    rightnode=leftnode
                    leftnode=max(rightnode-searchstep,Init_guess)
                    while(rightnode-leftnode>0.1):# once the feasible rightnode is found, the two separated method will be used to find the minimum feasbile solve
                        if ineqmycon((leftnode+rightnode)/2.0)==False:
                            leftnode=(leftnode+rightnode)/2.0
                        else:
                            rightnode=(leftnode+rightnode)/2.0
                else:
                    if leftnode>=rightnode:
                        if searchstep==0.01:
                            print("can not solve--Init_guess, y0,vy0,ay0,z0,vz0,az0,rightnode",Init_guess,y0,vy0,ay0,z0,vz0,az0,rightnode)
                            # rightnode=rightnode/increaseratio
                            leftnode=0.1
                            break
                        leftnode=Init_guess
                        searchstep=max(searchstep/2,0.01)
                    leftnode=leftnode+searchstep
                # print("searchstep--",searchstep)
                if time.time()-lastdroneupdatetime>=0.2:
                    print("Init_guess, y0,vy0,ay0,z0,vz0,az0",Init_guess,y0,vy0,ay0,z0,vz0,az0)
                    rightnode=rightnode/increaseratio # to avoid rightnode increase by multiply 1.5 many times in rightnode=rightnode*1.5
                    leftnode=0.1
                    print ("can not solve")
                    break

            running_time = time.time() - lastdroneupdatetime
            t=0
            if solveflag:
                lastsolved_time=time.time()
                print('time cost : %.5f sec' % running_time,'terminate cost : %.2f sec' % rightnode)
                # print('time cost : %.5f sec' % running_time,'leftnode cost : %.5f sec' % leftnode, 'right cost : %.5f sec' %rightnode,"y0,vy0,ay0,z0,vz0,az0",y0,vy0,ay0,z0,vz0,az0)
                # print('time cost : %.5f sec' % running_time,'terminate cost : %.5f sec' % Init_guess,"y0,vy0,ay0,z0,vz0,az0 ",y0,vy0,ay0,z0,vz0,az0)
                t=rightnode
                lastsolveduration=t
                # t=Init_guess
                if(t<=1):
                    increaseratio=1.0
            else:
                t=max(lastsolveduration-(time.time()-lastsolved_time),0)
                if(t>0.5):
                    decreasedtimemethod()
                    rightnode=t*increaseratio
                    leftnode=t*0.5
                    print('supplement decreased time----time cost : %.5f sec' % running_time,'terminate cost : %.2f sec' % t)
            if t>=0.5:
                controlstate_msg.discrepointpersecond = controlfreq
                controlstate_msg.arraylength = round(t*controlfreq)
                controlstate_msg.inicounter = (int)(min(2,controlstate_msg.arraylength))
                controlstate_msg.rendezvouswall_z=ztf_wall
                controlstate_msg.rendezvouswall_y=ytf_wall
                controlstate_msg.rendezvouswall_vy=vytf_wall
                controlstate_msg.rendezvouswall_vz=vztf_wall
                controlstate_msg.parabolictime=parabolictime


                times=np.linspace(0, 1, controlstate_msg.arraylength)*t
                tarray=np.array([[60/t**3,-360/t**4,720/t**5],[-24/t**2,168/t**3,-360/t**4],[3/t,-24/t**2,60/t**3]])

                alpha_y,beta_y,gamma_y=np.dot(tarray,np.array([aytf-ay0,vytf-vy0-ay0*t,ytf-y0-vy0*t-0.5*ay0*t**2]))
                alpha_z,beta_z,gamma_z=np.dot(tarray,np.array([aztf-az0,vztf-vz0-az0*t,ztf-z0-vz0*t-0.5*az0*t**2]))

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
                # print (y[-1],z[-1])
                currentupdateflag = False

                planned_path.header.stamp=rospy.Time.now()
                planned_path.header.frame_id="ground_link"
                planned_path.poses=[]

                for i in range(0, int(controlstate_msg.arraylength)):
                    planned_pose_stamped.pose.position.x=controlstate_msg.stateXarray[i]
                    planned_pose_stamped.pose.position.y=controlstate_msg.stateYarray[i]
                    planned_pose_stamped.pose.position.z=controlstate_msg.stateZarray[i]
                    planned_pose_stamped.header.stamp=rospy.Time.now()
                    # planned_pose_stamped.header.frame_id="base_link"
                    # planned_pose_stamped.header.seq=i
                    planned_path.poses.append(copy.deepcopy(planned_pose_stamped))

                path_pub.publish(planned_path)

        # print ("acc",az0)
        # recv_acc.linear_acceleration.y=ay0
        # recv_acc.linear_acceleration.z=az0
        # recv_acc.header.stamp=rospy.Time.now()
        # acc_pub.publish(recv_acc)

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
                # phiseries=-np.arctan2(ay,(az+9.8))
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
                # # print(result)
                # plt.show()

        rate.sleep()
        # print("sleep test")

if __name__ == '__main__':  # 主函数
    main()
