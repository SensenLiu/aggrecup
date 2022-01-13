#!/usr/bin/env python
# -*- coding: utf-8 -*-
# coding=utf-8
from __future__ import division
import copy
import socket
import numpy as np
from numpy import array
from numpy import dot
from numpy import sqrt
from numpy import diff
from numpy import arccos
from numpy import pi
from numpy import hstack
from numpy import append
from numpy import linspace
from numpy import full_like
from numpy import zeros_like
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
from geometry_msgs.msg import PointStamped
from geometry_msgs.msg import PoseStamped
from offb_posctl.msg import controlstate
from offb_posctl.msg import wallstate
from pyquaternion import Quaternion
# print(offb_posctl.__file__)

phi=1.22
normspeed=-0.5
tangentialspeed=0.3
xtf_wall=0.0
ytf_wall=10.0
ztf_wall=1.40 # true is 1.94
vxtf_wall=0.0
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
thrustmax=2*9.8
angleaccdmax=25
lbz=0.2
ubz=2.5
lbv=-5
ubv=5
currentupdateflag = False
wallstateupdateflag = False

wallstate_msg = wallstate()  # 要发布的控制量消息
lastwallupdatetime=0.0

lastsolved_time=0.0
lastsolveduration=0.0

waypoint=[[0.505578,0.619863,3.097399,-0.033222,-0.007129,0.805617],
          [2.088411,0.636662,3.532426,0.085445,3.882210,-0.092808],
          [4.166977,0.593732,4.446313,-0.154658,0.514369,0.832596],
          [6.406481,0.618986,4.573474,0.337406,1.566063,1.382526],
          [8.696024,0.961586,4.159587,1.029460,-3.553221,1.19400],
          [10,      1.4,     delta_vytf,delta_vztf,aytf, aztf]
]
waypoint=array(waypoint)
segmentstate=0
nearthreshold_y=3

def getterminateStateTime(t,solutionqualitycounter):
    global y0,xtf_wall,ytf_wall,ztf_wall,vxtf_wall,vytf_wall,vztf_wall,lastwallupdatetime, wallstate_msg,wallstateupdateflag
    global waypoint,segmentstate,nearthreshold_y
    temp_waypoint=waypoint.copy()

    if wallstateupdateflag:
        if((time.time()-lastwallupdatetime)*wallstate_msg.discrepointpersecond>=wallstate_msg.arraylength):
            wallstateupdateflag=False
            return False
        wallindex=int(min(round((t+time.time()-lastwallupdatetime+1.0/60.0)*wallstate_msg.discrepointpersecond),wallstate_msg.arraylength-1))#add 1/60.0 is to compliment the latency when offboard node receive the ocplan message,
        # this is because the offb_posctl_vicon node loops in 30Hz, and the receive will latency a period, therefore we add 1/60 as the middle of 1/30 to make sure the max time error located in 1/60
        # print(wallindex)
        if wallindex>=wallstate_msg.arraylength:
            return False

        xtf_wall=wallstate_msg.stateXarray[wallindex]
        vxtf_wall=wallstate_msg.stateVXarray[wallindex]
        ytf_wall=wallstate_msg.stateYarray[wallindex]
        vytf_wall=wallstate_msg.stateVYarray[wallindex]
        ztf_wall=wallstate_msg.stateZarray[wallindex]
        vztf_wall=wallstate_msg.stateVZarray[wallindex]


    temp_waypoint[:,0]=ytf_wall-waypoint[-1,0]+waypoint[:,0]-vytf_wall*t
    temp_waypoint[:,1]=ztf_wall-waypoint[-1,1]+waypoint[:,1]-vztf_wall*t
    temp_waypoint[:,2]=vytf_wall+waypoint[:,2]
    temp_waypoint[:,3]=vztf_wall+waypoint[:,3]
    temp_lbv=lbv+vytf_wall
    temp_ubv=ubv+vytf_wall
    pointcounter=1
    segmentstate=0
    while not (rospy.is_shutdown()) and pointcounter<=waypoint.shape[0]-1:
        if y0>=temp_waypoint[pointcounter-1,0]:
            segmentstate=pointcounter
            # print("segmentstate",segmentstate,'y0',y0)
        else:
            break
        pointcounter=pointcounter+1
    if segmentstate<=waypoint.shape[0]-2:
        if abs(temp_waypoint[segmentstate,0]-y0)<nearthreshold_y :
            ytf,ztf,vytf,vztf,aytf,aztf=temp_waypoint[segmentstate+1,:]
            t=t-(waypoint.shape[0]-2-segmentstate)*0.5
            # print("jump",'y0',y0,'temp_waypoint[segmentstate,0]',temp_waypoint[segmentstate,0] )
        else:
            ytf,ztf,vytf,vztf,aytf,aztf=temp_waypoint[segmentstate,:]
            t=t-(waypoint.shape[0]-1-segmentstate)*0.5
    else:
        ytf,ztf,vytf,vztf,aytf,aztf=temp_waypoint[segmentstate,:]
        t=t-(waypoint.shape[0]-1-segmentstate)*0.5

    # print("segment-----t----",t)

    return ytf,ztf,vytf,vztf,aytf,aztf,t,segmentstate,
# Constraint
def ineqmycon(x,solutionqualitycounter):
    global y0,z0,vy0,vz0,ay0,az0, ytf,ztf,vytf,vztf,aytf,aztf, thrustmax, angleaccdmax, lbz, ubz, lbv, ubv

    ytf,ztf,vytf,vztf,aytf,aztf,t,segmentstate=getterminateStateTime(x,solutionqualitycounter)
    # print('time in---time out',x,t,y0,z0,vy0,vz0,ay0,az0, ytf,ztf,vytf,vztf,aytf,aztf)
    if t<0:
        return False
    # t=1.57 # for debug
    # y0,z0,vy0,vz0,ay0,az0=[-0.06390421092510223, 0.53312087059021, 0.028850866481661797, -0.07197022438049316, 0.021810559555888176, 0.010911194607615471]
    # ytf,ztf,vytf,vztf,aytf,aztf=[0.505578, 0.619863, 3.097399, -0.033222, -0.007129, 0.805617]
    # print("t----wallindex----ytf:",t,wallindex,ytf)
    tarray=array([[60/t**3,-360/t**4,720/t**5],[-24/t**2,168/t**3,-360/t**4],[3/t,-24/t**2,60/t**3]])
    alpha_y,beta_y,gamma_y=dot(tarray,array([aytf-ay0,vytf-vy0-ay0*t,ytf-y0-vy0*t-0.5*ay0*t**2]))
    alpha_z,beta_z,gamma_z=dot(tarray,array([aztf-az0,vztf-vz0-az0*t,ztf-z0-vz0*t-0.5*az0*t**2]))
    tmesh=t*(array(linspace(0.001, 1, 20)))

    ## z's lower bound  constraints
    y=alpha_y/120*tmesh**5+beta_y/24*tmesh**4+gamma_y/6*tmesh**3+ay0/2*tmesh**2+vy0*tmesh+y0
    z=alpha_z/120*tmesh**5+beta_z/24*tmesh**4+gamma_z/6*tmesh**3+az0/2*tmesh**2+vz0*tmesh+z0
    c2=lbz-z
    if(sum(c2>0.05)>0):
        return False

    ##to judge the curve is located on the identical side of the straight line (y0,z0)--(ytf,ztf)
    c2_2=(y-y0)*(ytf-y0)+(z-z0)*(ztf-z0)
    dotnormal=sqrt(((y-y0)**2+(z-z0)**2)*((ztf-z0)**2+(ytf-y0)**2))
    # c2_2=np.delete(c2_2,np.where(dotnormal==0))
    # dotnormal=np.delete(dotnormal,np.where(dotnormal==0))
    c2_3=c2_2/dotnormal
    c2_3=np.delete(c2_3,np.where(c2_3>1))
    c2_3=np.delete(c2_3,np.where(c2_3<-1))
    # print("c2_3",c2_3)
    c2_3=diff(arccos(c2_3)*180/pi)
    if(sum(c2_3>0.05)>0):# the first one or three points at the begnining usually are outliers and it can be accepted
        return False

    ##velocity constraints
    vy=alpha_y/24*tmesh**4+beta_y/6*tmesh**3+gamma_y/2*tmesh**2+ay0*tmesh+vy0
    vz=alpha_z/24*tmesh**4+beta_z/6*tmesh**3+gamma_z/2*tmesh**2+az0*tmesh+vz0
    c10=vy-ubv
    c11=lbv-vy
    c12=vz-ubv
    c13=lbv-vz
    if(sum(hstack((c10,c11,c12,c13))>0.05)>0):
        return False

    ##actuator constraints
    ay=alpha_y/6*tmesh**3+beta_y/2*tmesh**2+gamma_y*tmesh+ay0
    az=alpha_z/6*tmesh**3+beta_z/2*tmesh**2+gamma_z*tmesh+az0
    jy=alpha_y/2*tmesh**2+beta_y*tmesh+gamma_y
    jz=alpha_z/2*tmesh**2+beta_z*tmesh+gamma_z
    numerator=(jy*(az+9.8)-ay*jz)
    denominator=((az+9.8)**2+ay**2)
    numerator=np.delete(numerator,np.where(denominator==0))
    denominator=np.delete(denominator,np.where(denominator==0))
    angleacc=append(0,-diff(numerator/denominator)/(tmesh[2]-tmesh[1]))

    # angleacc=((((alpha_y*tmesh**2)/2 + beta_y*tmesh + gamma_y)/((alpha_z*tmesh**3)/6 + (beta_z*tmesh**2)/2 + gamma_z*tmesh + az0 + 49/5) -
    #            (((alpha_z*tmesh**2)/2 + beta_z*tmesh + gamma_z)*((alpha_y*tmesh**3)/6 + (beta_y*tmesh**2)/2 + gamma_y*tmesh + ay0))/
    #            ((alpha_z*tmesh**3)/6 + (beta_z*tmesh**2)/2 + gamma_z*tmesh + az0 + 49/5)**2)*((2*((alpha_y*tmesh**2)/2 + beta_y*tmesh + gamma_y)*
    #            ((alpha_y*tmesh**3)/6 + (beta_y*tmesh**2)/2 + gamma_y*tmesh + ay0))/((alpha_z*tmesh**3)/6 + (beta_z*tmesh**2)/2 + gamma_z*tmesh +
    #            az0 + 49/5)**2 - (2*((alpha_z*tmesh**2)/2 + beta_z*tmesh + gamma_z)*((alpha_y*tmesh**3)/6 + (beta_y*tmesh**2)/2 + gamma_y*tmesh + ay0)**2)/
    #            ((alpha_z*tmesh**3)/6 + (beta_z*tmesh**2)/2 + gamma_z*tmesh + az0 + 49/5)**3))/(((alpha_y*tmesh**3)/6 + (beta_y*tmesh**2)/2 + gamma_y*tmesh + ay0)**2/
    #            ((alpha_z*tmesh**3)/6 + (beta_z*tmesh**2)/2 + gamma_z*tmesh + az0 + 49/5)**2 + 1)**2 - ((beta_y + alpha_y*tmesh)/((alpha_z*tmesh**3)/6 +
    #            (beta_z*tmesh**2)/2 + gamma_z*tmesh + az0 + 49/5) - ((beta_z + alpha_z*tmesh)*((alpha_y*tmesh**3)/6 + (beta_y*tmesh**2)/2 + gamma_y*tmesh +
    #            ay0))/((alpha_z*tmesh**3)/6 + (beta_z*tmesh**2)/2 + gamma_z*tmesh + az0 + 49/5)**2 - (2*((alpha_y*tmesh**2)/2 + beta_y*tmesh + gamma_y)*
    #            ((alpha_z*tmesh**2)/2 + beta_z*tmesh + gamma_z))/((alpha_z*tmesh**3)/6 + (beta_z*tmesh**2)/2 + gamma_z*tmesh + az0 + 49/5)**2 +
    #            (2*((alpha_z*tmesh**2)/2 + beta_z*tmesh + gamma_z)**2*((alpha_y*tmesh**3)/6 + (beta_y*tmesh**2)/2 + gamma_y*tmesh + ay0))/((alpha_z*tmesh**3)/6 +
    #            (beta_z*tmesh**2)/2 + gamma_z*tmesh + az0 + 49/5)**3)/(((alpha_y*tmesh**3)/6 + (beta_y*tmesh**2)/2 + gamma_y*tmesh + ay0)**2/
    #            ((alpha_z*tmesh**3)/6 + (beta_z*tmesh**2)/2 + gamma_z*tmesh + az0 + 49/5)**2 + 1)

    thrust=sqrt((ay**2 + (az + 49/5)**2))
    c3=-angleacc*thrustmax/(4*angleaccdmax)+thrust/2-9.8
    c4=angleacc*thrustmax/(4*angleaccdmax)-thrust/2
    c5=angleacc*thrustmax/(4*angleaccdmax)+thrust/2-9.8
    c6=-angleacc*thrustmax/(4*angleaccdmax)-thrust/2
    if(sum(hstack((c3,c4,c5,c6))>0.05)>0):
        return False

    return True

def pos_twist_callback(data):
    global vy0, y0, vz0, z0, ay0,az0,currentupdateflag,targetstartmoveflag
    y0,z0,vy0,vz0,ay0,az0 = data.pose.pose.position.y,data.pose.pose.position.z,data.twist.twist.linear.y,data.twist.twist.linear.z,data.twist.twist.angular.y,data.twist.twist.angular.z
    # print(y0)
    currentupdateflag = True


def droneImu_callback(data):
    global ay0, az0
    # ay0 = data.linear_acceleration.y
    # az0 = data.linear_acceleration.z-9.8
    # qcurrent = Quaternion(data.orientation.w, data.orientation.x, data.orientation.y, data.orientation.z)
    # qaccinworld=qcurrent.rotate(Quaternion(vector=[data.linear_acceleration.x,data.linear_acceleration.y,data.linear_acceleration.z]))
    # ay0=0.2*min(max(qaccinworld.y,-2*9.8),2*9.8)+0.8*ay0
    # az0=0.2*min(max(qaccinworld.z-9.8,-9.8),9.8)+0.8*az0
    # ay0=0
    # az0=0
    # print ("ay0",ay0)

def wallstate_callback(data):
    global wallstate_msg,wallstateupdateflag,lastwallupdatetime
    wallstate_msg=data
    wallstateupdateflag=True
    # print("time.time(): ",time.time())
    # print("wallstate_msg.header.stamp--100---: ",wallstate_msg.header.stamp.secs+wallstate_msg.header.stamp.nsecs/1e9)
    lastwallupdatetime=wallstate_msg.header.stamp.secs+wallstate_msg.header.stamp.nsecs/1e9

def main():
    global currentupdateflag, wallstateupdateflag, targetstartmoveflag,lastsolved_time,lastsolveduration,xtf_wall,ytf_wall,ztf_wall,vxtf_wall,vytf_wall,vztf_wall,normspeed,wallstate_msg
    controlstate_msg = controlstate()  # 要发布的控制量消息
    planned_path=Path()
    target_path=Path()
    planned_pose_stamped=PoseStamped()
    target_pose_stamped=PoseStamped()
    waypoints=PointStamped()
    recv_acc=Imu()

    rospy.init_node('minimumsnap_control', anonymous=True)
    uav_id = rospy.get_param("~id", "")
    rate = rospy.Rate(100)

    rospy.Subscriber(uav_id + "current_relative_postwist",
                     Odometry, pos_twist_callback,buff_size=2**24,tcp_nodelay=True)
    # rospy.Subscriber(uav_id + "/mavros/imu/data",
    #                  Imu, droneImu_callback)
    rospy.Subscriber(uav_id + "predictedwall_array",
                     wallstate, wallstate_callback, queue_size=1,buff_size=2**24,tcp_nodelay=True)  # plane velocity
    control_pub = rospy.Publisher(uav_id + "ocp/control_state", controlstate, queue_size=1)
    path_pub = rospy.Publisher(uav_id + "plannedtrajectory",Path,queue_size=1)
    point_pub = rospy.Publisher(uav_id + "waypoints",PointStamped,queue_size=1)
    acc_pub = rospy.Publisher(uav_id +"/acc_data_inworld",Imu,queue_size=1)

    # currentupdateflag=True # if we use True as initial value, the initial terminate time will be less than the firt true time, and lead no resolution always
    currentupdateflag=False
    ytf, vytf,aytf,ztf,vztf,aztf=0,0,0,0,0,0
    y,vy,ay,z,vz,az,lastarraylength=0,0,0,0,0,0,0
    controlfreq=30
    planstoptime=0.2
    Init_guess,increaseratio=0.1,1.5
    searchstep,leftnode,rightnode,solutionqualitycounter,segmentstate,solveflag=1, 0.1, 5, 0,0,False
    rospy.loginfo_throttle(0.02,("Terminate state, ytf,vytf,aytf,ztf,vztf,aztf : %.2f  %.2f %.2f %.2f %.2f %.2f",ytf,vytf,aytf,ztf,vztf,aztf))

    # while not((rospy.is_shutdown())): ## test the runtime of the core function
    #     startsolvetime=time.time()
    #     ineqmycon(10,solutionqualitycounter)
    #     running_time = time.time() - startsolvetime
    #     print(running_time)

    while not((rospy.is_shutdown())) and wallstateupdateflag:
        if wallstateupdateflag:
            if math.fabs(wallstate_msg.stateVYarray[0])>0.1:
                break
        rospy.loginfo_throttle(0.02,"the wall have not been moved_vytf_wall:")
        if wallstateupdateflag:
            rospy.loginfo_throttle(0.02,wallstate_msg.stateVYarray[0])
        rate.sleep()

    print("phi,cuplength~~~~~",phi,cuplength)
    for j in range(0,10):## display the waypoint. The "History Length" should be not less than 5 in rivz pointstameped property configurateion
        for i in range(0, 6):
            waypoints.header.stamp=rospy.Time.now()
            waypoints.header.frame_id="ground_link"
            waypoints.point.x=0
            waypoints.point.y=waypoint[i,0]
            waypoints.point.z=waypoint[i,1]
            # print(waypoints.point.y)
            point_pub.publish(waypoints)
            rate.sleep()

    while not (rospy.is_shutdown()):
        # if currentupdateflag and wallstateupdateflag:
        if currentupdateflag:
            Init_guess=leftnode*0.5
            leftnode=Init_guess
            rightnode=min(rightnode*increaseratio,5)
            searchstep=max((rightnode-leftnode)/5, 0.01)
            solveflag=False
            startsolvetime=time.time()

            while rightnode-leftnode>0.2 or solveflag==False:
                solveflag=ineqmycon(leftnode,solutionqualitycounter)
                # print("solveflag,leftnode",solveflag,leftnode)
                if solveflag==True:
                    rightnode=leftnode
                    leftnode=max(rightnode-searchstep,Init_guess)
                    while(rightnode-leftnode>0.1):# once the feasible rightnode is found, the two separated method will be used to find the minimum feasbile solve
                        if ineqmycon((leftnode+rightnode)/2.0,solutionqualitycounter)==False:
                            leftnode=(leftnode+rightnode)/2.0
                        else:
                            rightnode=(leftnode+rightnode)/2.0
                    break # this break is necessary, becasue it is possible that "if time.time()-startsolvetime>=0.03" will work,
                    # the after the while and this make the sentence " rightnode=rightnode/increaseratio" executed and reduce the rightnode
                    # and then make the t is negative in "if soveflag:"part
                else:
                    if leftnode>=rightnode:
                        if searchstep==0.05:
                            ytf,ztf,vytf,vztf,aytf,aztf,t,segmentstate=getterminateStateTime(rightnode,solutionqualitycounter)
                            # print('cannot solve-Init_guess %.2f' % Init_guess,  'y0:%.2f' % y0,  'vy0:%.2f' %vy0,  'ay0:%.2f' % ay0,   'z0:%.2f' % z0,   'vz0:%.2f' % vz0,
                            #     'az0:%.2f' %az0,  'ytf:%.2f' % ytf,  'vytf:%.2f' % vytf,  'aytf:%.2f' % aytf,   'ztf:%.2f' % ztf,   'vztf:%.2f' % vztf, 'aztf:%.2f' % aztf,
                            #     'rightnode %.2f' % rightnode, 't:%.2f' % t, 'segmentstate %.1f' % segmentstate)
                            rospy.loginfo_throttle(0.02,('cannot solve-Init_guess %.2f' % Init_guess,  'y0: %.2f' % y0, 'z0:%.2f' % z0,
                            'ytf:%.2f' % ytf,  'vyft: %.2f' % vytf, 'ztf:%.2f' % ztf, 'rightnode %.2f' % rightnode, 't:%.2f' % t, 'segmentstate %d' % segmentstate))
                            leftnode=0.1
                            break
                        leftnode=Init_guess
                        searchstep=max(searchstep/2,0.01)
                    leftnode=min(leftnode+searchstep,rightnode)
                # print("searchstep--",searchstep)
                if time.time()-startsolvetime>=0.03:
                    ytf,ztf,vytf,vztf,aytf,aztf,t,segmentstate=getterminateStateTime(rightnode,solutionqualitycounter)
                    # print('solvetimeout-Init_guess %.2f' % Init_guess,  'y0:%.2f' % y0,  'vy0:%.2f' %vy0,  'ay0:%.2f' % ay0,   'z0:%.2f' % z0,   'vz0:%.2f' % vz0,
                    #       'az0:%.2f' %az0,  'ytf:%.2f' % ytf,  'vytf:%.2f' % vytf,  'aytf:%.2f' % aytf,   'ztf:%.2f' % ztf,   'vztf:%.2f' % vztf, 'aztf:%.2f' % aztf,
                    #       't:%.2f' % t, 'segmentstate %.1f' % segmentstate)
                    rospy.loginfo_throttle(0.02,('timeout solve-Init_guess %.2f' % Init_guess,  'y0:%.2f' % y0, 'z0:%.2f' % z0,
                     'ytf:%.2f' % ytf,  'ztf:%.2f' % ztf, 'vyft: %.2f' % vytf, 'rightnode %.2f' % rightnode, 't:%.2f' % t, 'segmentstate %d' % segmentstate))
                    rightnode=rightnode/increaseratio # to avoid rightnode increase by multiply 1.5 many times in rightnode=rightnode*1.5
                    leftnode=0.1
                    break

            running_time = time.time() - startsolvetime
            segmenttime=0
            if solveflag:
                lastsolved_time=time.time()
                lastsolveduration=rightnode-(controlstate_msg.inicounter/controlfreq)
                residualduration=lastsolveduration
                ytf,ztf,vytf,vztf,aytf,aztf,t,segmentstate=getterminateStateTime(lastsolveduration,solutionqualitycounter)# this function is very important. because when rightnode-leftnode<0.1 and the loop break.
                # it is very probabily that the t is left node and ineqmycon(leftnode) will change the teriminate condition to former time point.
                # Therefore, the getterminateStateTime() fucntion use the rightnode time to calculate the rendezvous point again and recover the
                # terminal condition represented by global arguments in ineqmycon() fucntion
                # print('controlstate_msg.inicounter/controlfreq : %.5f sec' % (controlstate_msg.inicounter/controlfreq))
                segmenttime=t+(controlstate_msg.inicounter/controlfreq)# recovery t to the duration corresponding planning start state
                rospy.loginfo_throttle(0.02,('solved-computercost : %.5f' % running_time, 'segmenttime : %.2f' % segmenttime, 'segment %d' % segmentstate,
                'residualduration %.2f' % residualduration,'y0: %.3f' % y0,'ytf: %.3f' % ytf, 'ztf : %.2f' % ztf,'vyft: %.2f' % vytf))
                solutionqualitycounter=0
            else:
                residualduration=max(lastsolveduration-(time.time()-lastsolved_time)-(controlstate_msg.inicounter/controlfreq),0)
                if(residualduration>planstoptime):
                    ytf,ztf,vytf,vztf,aytf,aztf,t,segmentstate=getterminateStateTime(residualduration,solutionqualitycounter)
                    if t>0:
                        segmenttime=t+(controlstate_msg.inicounter/controlfreq)# recovery t to the duration corresponding planning start state
                        leftnode,rightnode=residualduration,residualduration ## for next solve
                        rospy.loginfo_throttle(0.02,('supple-computercost : %.5f' % running_time, 'segmenttime : %.2f' % segmenttime, 'segment %d' % segmentstate,
                        'residualduration %.2f' % residualduration,'y0: %.3f' % y0,'ytf: %.3f' % ytf, 'ztf : %.2f' % ztf,'vyft: %.2f' % vytf))
                solutionqualitycounter=min(solutionqualitycounter+1,3)
                print('solutionqualitycounter',solutionqualitycounter)
            if(residualduration<=1.1 and residualduration>0):
                increaseratio=1.0

            if residualduration>planstoptime:
                t=segmenttime
                if solveflag: ## solved success
                    controlstate_msg.header.stamp=rospy.Time.now()
                    controlstate_msg.discrepointpersecond = controlfreq
                    controlstate_msg.arraylength = round(t*controlfreq)+1# this "add 1" is very important, because if we want n intervals,we should have n+1 array length.
                    controlstate_msg.inicounter = (int)(min(2,controlstate_msg.arraylength))
                    controlstate_msg.timecost=running_time
                    controlstate_msg.increaseratio=increaseratio
                    controlstate_msg.segment=segmentstate
                    lastarraylength=controlstate_msg.arraylength
                    if wallstateupdateflag:
                        controlstate_msg.currentwall_x=wallstate_msg.currentwall_x
                        controlstate_msg.currentwall_y=wallstate_msg.currentwall_y
                        controlstate_msg.currentwall_z=wallstate_msg.currentwall_z
                        controlstate_msg.currentwall_vx=wallstate_msg.currentwall_vx
                        controlstate_msg.currentwall_vy=wallstate_msg.currentwall_vy
                        controlstate_msg.currentwall_vz=wallstate_msg.currentwall_vz

                        controlstate_msg.predictstartwall_x=wallstate_msg.stateXarray[0]
                        controlstate_msg.predictstartwall_y=wallstate_msg.stateYarray[0]
                        controlstate_msg.predictstartwall_z=wallstate_msg.stateZarray[0]
                        controlstate_msg.predictstartwall_vx=wallstate_msg.stateVXarray[0]
                        controlstate_msg.predictstartwall_vy=wallstate_msg.stateVYarray[0]
                        controlstate_msg.predictstartwall_vz=wallstate_msg.stateVZarray[0]
                        # print("wallstate_msg.header.stamp---: ",wallstate_msg.header.stamp.secs+wallstate_msg.header.stamp.nsecs/1e9)
                        # print('wallstate_msg.currentwall_y : %.5f sec' % wallstate_msg.currentwall_y,'wallstate_msg.stateYarray[0] : %.5f sec' % wallstate_msg.stateYarray[0])
                    else:
                        controlstate_msg.currentwall_x=xtf_wall
                        controlstate_msg.currentwall_y=ytf_wall
                        controlstate_msg.currentwall_z=ztf_wall
                        controlstate_msg.currentwall_vx=vxtf_wall
                        controlstate_msg.currentwall_vy=vytf_wall
                        controlstate_msg.currentwall_vz=vztf_wall
                        controlstate_msg.predictstartwall_x=xtf_wall
                        controlstate_msg.predictstartwall_y=ytf_wall
                        controlstate_msg.predictstartwall_z=ztf_wall
                        controlstate_msg.predictstartwall_vx=vxtf_wall
                        controlstate_msg.predictstartwall_vy=vytf_wall
                        controlstate_msg.predictstartwall_vz=vztf_wall

                    controlstate_msg.rendezvouswall_x=xtf_wall
                    controlstate_msg.rendezvouswall_y=ytf_wall
                    controlstate_msg.rendezvouswall_z=ztf_wall
                    controlstate_msg.rendezvouswall_vx=vxtf_wall
                    controlstate_msg.rendezvouswall_vy=vytf_wall
                    controlstate_msg.rendezvouswall_vz=vztf_wall
                    controlstate_msg.parabolictime=parabolictime
                    # print("controlstate_msg.rendezvouswall_y", controlstate_msg.rendezvouswall_y)

                    # print("Terminate state before calcualted, t,ytf,vytf,aytf,ztf,vztf,aztf,y0,vy0,ay0,z0,vz0,az0",t,ytf,vytf,aytf,ztf,vztf,aztf,y0,vy0,ay0,z0,vz0,az0)

                    times=linspace(0, 1, controlstate_msg.arraylength)*t
                    tarray=array([[60/t**3,-360/t**4,720/t**5],[-24/t**2,168/t**3,-360/t**4],[3/t,-24/t**2,60/t**3]])

                    alpha_y,beta_y,gamma_y=dot(tarray,array([aytf-ay0,vytf-vy0-ay0*t,ytf-y0-vy0*t-0.5*ay0*t**2]))
                    alpha_z,beta_z,gamma_z=dot(tarray,array([aztf-az0,vztf-vz0-az0*t,ztf-z0-vz0*t-0.5*az0*t**2]))

                    y=alpha_y/120*times**5+beta_y/24*times**4+gamma_y/6*times**3+ay0/2*times**2+vy0*times+y0
                    vy=alpha_y/24*times**4+beta_y/6*times**3+gamma_y/2*times**2+ay0*times+vy0
                    ay=alpha_y/6*times**3+beta_y/2*times**2+gamma_y*times+ay0

                    z=alpha_z/120*times**5+beta_z/24*times**4+gamma_z/6*times**3+az0/2*times**2+vz0*times+z0
                    vz=alpha_z/24*times**4+beta_z/6*times**3+gamma_z/2*times**2+az0*times+vz0
                    az=alpha_z/6*times**3+beta_z/2*times**2+gamma_z*times+az0

                    controlstate_msg.stateXarray = full_like(times,controlstate_msg.rendezvouswall_x)
                    controlstate_msg.stateYarray = y
                    controlstate_msg.stateZarray = z
                    controlstate_msg.stateVXarray = full_like(times,controlstate_msg.rendezvouswall_vx)
                    controlstate_msg.stateVYarray = vy
                    controlstate_msg.stateVZarray = vz
                    controlstate_msg.stateAXarray = zeros_like(times)
                    controlstate_msg.stateAYarray = ay
                    controlstate_msg.stateAZarray = az
                    # print("y--------", y)
                else: ## cannot solve
                    controlstate_msg.header.stamp=rospy.Time.now()
                    controlstate_msg.discrepointpersecond = controlfreq
                    controlstate_msg.arraylength = round(t*controlfreq)+1# this "add 1" is very important, because if we want n intervals,we should have n+1 array length.
                    controlstate_msg.inicounter = (int)(min(2,controlstate_msg.arraylength))
                    controlstate_msg.timecost=running_time
                    controlstate_msg.increaseratio=increaseratio
                    controlstate_msg.segment=segmentstate

                    if wallstateupdateflag:
                        controlstate_msg.currentwall_x=wallstate_msg.currentwall_x
                        controlstate_msg.currentwall_y=wallstate_msg.currentwall_y
                        controlstate_msg.currentwall_z=wallstate_msg.currentwall_z
                        controlstate_msg.currentwall_vx=wallstate_msg.currentwall_vx
                        controlstate_msg.currentwall_vy=wallstate_msg.currentwall_vy
                        controlstate_msg.currentwall_vz=wallstate_msg.currentwall_vz

                        controlstate_msg.predictstartwall_x=wallstate_msg.stateXarray[0]
                        controlstate_msg.predictstartwall_y=wallstate_msg.stateYarray[0]
                        controlstate_msg.predictstartwall_z=wallstate_msg.stateZarray[0]
                        controlstate_msg.predictstartwall_vx=wallstate_msg.stateVXarray[0]
                        controlstate_msg.predictstartwall_vy=wallstate_msg.stateVYarray[0]
                        controlstate_msg.predictstartwall_vz=wallstate_msg.stateVZarray[0]
                        # print("wallstate_msg.header.stamp---: ",wallstate_msg.header.stamp.secs+wallstate_msg.header.stamp.nsecs/1e9)
                        # print('wallstate_msg.currentwall_y : %.5f sec' % wallstate_msg.currentwall_y,'wallstate_msg.stateYarray[0] : %.5f sec' % wallstate_msg.stateYarray[0])
                    else:
                        controlstate_msg.currentwall_x=xtf_wall
                        controlstate_msg.currentwall_y=ytf_wall
                        controlstate_msg.currentwall_z=ztf_wall
                        controlstate_msg.currentwall_vx=vxtf_wall
                        controlstate_msg.currentwall_vy=vytf_wall
                        controlstate_msg.currentwall_vz=vztf_wall
                        controlstate_msg.predictstartwall_x=xtf_wall
                        controlstate_msg.predictstartwall_y=ytf_wall
                        controlstate_msg.predictstartwall_z=ztf_wall
                        controlstate_msg.predictstartwall_vx=vxtf_wall
                        controlstate_msg.predictstartwall_vy=vytf_wall
                        controlstate_msg.predictstartwall_vz=vztf_wall

                    controlstate_msg.rendezvouswall_x=xtf_wall
                    controlstate_msg.rendezvouswall_y=ytf_wall
                    controlstate_msg.rendezvouswall_z=ztf_wall
                    controlstate_msg.rendezvouswall_vx=vxtf_wall
                    controlstate_msg.rendezvouswall_vy=vytf_wall
                    controlstate_msg.rendezvouswall_vz=vztf_wall
                    controlstate_msg.parabolictime=parabolictime
                    interceptindex=int(max(min(lastarraylength-round(t*controlfreq)-1,lastarraylength-1),0))

                    controlstate_msg.stateXarray = full_like(y[interceptindex:],controlstate_msg.rendezvouswall_x)
                    controlstate_msg.stateYarray = y[interceptindex:]+ytf-y[-1]
                    controlstate_msg.stateZarray = z[interceptindex:]+ztf-z[-1]
                    controlstate_msg.stateVXarray = full_like(vy[interceptindex:],controlstate_msg.rendezvouswall_vx)
                    controlstate_msg.stateVYarray = vy[interceptindex:]
                    controlstate_msg.stateVZarray = vz[interceptindex:]
                    controlstate_msg.stateAXarray = zeros_like(ay[interceptindex:])
                    controlstate_msg.stateAYarray = ay[interceptindex:]
                    controlstate_msg.stateAZarray = az[interceptindex:]

                control_pub.publish(controlstate_msg)
                # print (y[-1],z[-1])

                planned_path.header.stamp=rospy.Time.now()
                planned_path.header.frame_id="ground_link"
                planned_path.poses=[]

                for i in range(0, int(len(controlstate_msg.stateXarray))):
                    planned_pose_stamped.pose.position.x=controlstate_msg.stateXarray[i]
                    planned_pose_stamped.pose.position.y=controlstate_msg.stateYarray[i]
                    planned_pose_stamped.pose.position.z=controlstate_msg.stateZarray[i]
                    planned_pose_stamped.header.stamp=rospy.Time.now()
                    # planned_pose_stamped.header.frame_id="base_link"
                    # planned_pose_stamped.header.seq=i
                    planned_path.poses.append(copy.deepcopy(planned_pose_stamped))

                path_pub.publish(planned_path)
            currentupdateflag = False

        recv_acc.linear_acceleration.y=ay0
        recv_acc.linear_acceleration.z=az0
        recv_acc.header.stamp=rospy.Time.now()
        acc_pub.publish(recv_acc)
        # rate.sleep()#do not use it to make sure the plan can be run immediately

if __name__ == '__main__':  # 主函数
    main()
