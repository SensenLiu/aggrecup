#!/usr/bin/env python
# -*- coding: utf-8 -*-
# coding=utf-8
import socket
import numpy as np
from scipy.optimize import minimize
import time
from numba import jit, float64
import datetime
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy
print(np.__version__)

phi=1.22
normspeed=-0.5
tangentialspeed=0.3
ytf_wall=2.00
ztf_wall=1.40 # true is 1.94
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

lastwallupdatetime=0.0
lastdroneupdatetime=0.0

lastsolved_time=0.0
lastsolveduration=0.0

y0,z0,vy0,vz0,ay0,az0=[8.119218826293945, 1.5682785511016846, 5.646708965301514, 0.9292353391647339, 6.371948719024658, -5.756661415100098]
ytf,ztf,vytf,vztf,aytf,aztf=[10.867484919410352, 1.2385329108834278, 6.950641729776516, -0.15108810566812603, 0.514369, 0.832596]
print("Terminate state, ytf,vytf,aytf,ztf,vztf,aztf",ytf,vytf,aytf,ztf,vztf,aztf)

def ineqmycon(x):
    global ay0, vy0, y0, az0, vz0, z0, aytf, vytf, ytf, aztf, vztf, ztf, meshpoint, thrustmax, angleaccdmax, lbz, lbv, ubv
    t=x[0]
    tarray=np.array([[60/t**3,-360/t**4,720/t**5],[-24/t**2,168/t**3,-360/t**4],[3/t,-24/t**2,60/t**3]])

    alpha_y,beta_y,gamma_y=np.dot(tarray,np.array([aytf-ay0,vytf-vy0-ay0*t,ytf-y0-vy0*t-0.5*ay0*t**2]))
    alpha_z,beta_z,gamma_z=np.dot(tarray,np.array([aztf-az0,vztf-vz0-az0*t,ztf-z0-vz0*t-0.5*az0*t**2]))


    tmesh=t*(np.array(meshpoint))
    angleacc=np.zeros_like(tmesh)
    for i in range(len(tmesh)):
        # print("i===",tmesh)
        t=tmesh[i]
        angleacc[i]=((((alpha_y*t**2)/2 + beta_y*t + gamma_y)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - (((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2)*((2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3))/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)**2 - ((beta_y + alpha_y*t)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - ((beta_z + alpha_z*t)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_z*t**2)/2 + beta_z*t + gamma_z))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)**2*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3)/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)

    t=x[0]
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
    # c7=((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 +9.8) # it is not necessary, because the specified terminal acc will locate in the feasible region
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
    # tmesh=t*(np.array(meshpoint))
    # angleacc=np.zeros_like(tmesh)
    # for i in range(len(tmesh)):
    #     # print("i===",tmesh)
    #     t=tmesh[i]
    #     angleacc[i]=((((alpha_y*t**2)/2 + beta_y*t + gamma_y)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - (((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2)*((2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3))/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)**2 - ((beta_y + alpha_y*t)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - ((beta_z + alpha_z*t)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_z*t**2)/2 + beta_z*t + gamma_z))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)**2*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3)/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)


    c10=-(alpha_y/24*tmesh**4+beta_y/6*tmesh**3+gamma_y/2*tmesh**2+ay0*tmesh+vy0-ubv)
    c11=-(lbv-(alpha_y/24*tmesh**4+beta_y/6*tmesh**3+gamma_y/2*tmesh**2+ay0*tmesh+vy0))
    c12=-(alpha_z/24*tmesh**4+beta_z/6*tmesh**3+gamma_z/2*tmesh**2+az0*tmesh+vz0-ubv)
    c13=-(lbv-(alpha_z/24*tmesh**4+beta_z/6*tmesh**3+gamma_z/2*tmesh**2+az0*tmesh+vz0))
    # print("--------", t,np.hstack((c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13)))
    # print("--------", (np.hstack((c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13))>0).all())
    return (np.hstack((c2,c3,c4,c5,c6,c10,c11,c12,c13))>-0.05).all()

def main():
    startsolvetime=time.time()
    ineqmycon(np.array([1]))
    running_time = time.time() - startsolvetime
    print(running_time)
    Initial_guess=np.array([0.1])
    # while True:
    Initial_guess=Initial_guess*0.5
    start = time.time()
    start0=time.time()
    resetcounter=0
    controlfreq=30
    while (ineqmycon(Initial_guess)==False):
        Initial_guess=Initial_guess+0.01
        # print(ineqmycon(Initial_guess))
        end = time.time()
        if((end-start)>0.5):
            if(resetcounter==0):
                Initial_guess=np.array([0.1])
                start = time.time()
                resetcounter=1
            else:
                print("can not solve")
                break

    running_time = end - start0
    print('time cost : %.5f sec' % running_time)
    print(Initial_guess[0], ineqmycon(Initial_guess))
    Initial_guess[0]=0.4
    times=np.linspace(0,1,round(Initial_guess[0]*controlfreq))*Initial_guess

    t=Initial_guess[0]
    y0,z0,vy0,vz0,ay0,az0=[8.119218826293945, 1.5682785511016846, 5.646708965301514, 0.9292353391647339, 6.371948719024658, -5.756661415100098]
    ytf,ztf,vytf,vztf,aytf,aztf=[10.867484919410352, 1.2385329108834278, 6.950641729776516, -0.15108810566812603, 0.514369, 0.832596]
    tarray=np.array([[60/t**3,-360/t**4,720/t**5],[-24/t**2,168/t**3,-360/t**4],[3/t,-24/t**2,60/t**3]])

    tarray=np.array([[60/t**3,-360/t**4,720/t**5],[-24/t**2,168/t**3,-360/t**4],[3/t,-24/t**2,60/t**3]])

    alpha_y,beta_y,gamma_y=np.dot(tarray,np.array([aytf-ay0,vytf-vy0-ay0*t,ytf-y0-vy0*t-0.5*ay0*t**2]))
    alpha_z,beta_z,gamma_z=np.dot(tarray,np.array([aztf-az0,vztf-vz0-az0*t,ztf-z0-vz0*t-0.5*az0*t**2]))

    y=alpha_y/120*times**5+beta_y/24*times**4+gamma_y/6*times**3+ay0/2*times**2+vy0*times+y0
    vy=alpha_y/24*times**4+beta_y/6*times**3+gamma_y/2*times**2+ay0*times+vy0
    ay=alpha_y/6*times**3+beta_y/2*times**2+gamma_y*times+ay0

    z=alpha_z/120*times**5+beta_z/24*times**4+gamma_z/6*times**3+az0/2*times**2+vz0*times+z0
    vz=alpha_z/24*times**4+beta_z/6*times**3+gamma_z/2*times**2+az0*times+vz0
    az=alpha_z/6*times**3+beta_z/2*times**2+gamma_z*times+az0

    a=np.sqrt(az**2+ay**2)
    thurst=np.sqrt((az+9.8)**2+ay**2)
    phiseries=-np.arctan2(ay,(az+9.8))
    # print("az--------", az)
    angleacc=np.zeros_like(times)
    for i in range(len(times)):
        t=times[i]
        angleacc[i]=((((alpha_y*t**2)/2 + beta_y*t + gamma_y)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - (((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2)*((2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3))/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)**2 - ((beta_y + alpha_y*t)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - ((beta_z + alpha_z*t)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_z*t**2)/2 + beta_z*t + gamma_z))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)**2*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3)/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)

    F1=angleacc*thrustmax/(4*angleaccdmax)+thurst/2
    F2=-angleacc*thrustmax/(4*angleaccdmax)+thurst/2

    plotlinewidth = 2
    plotfontsize = 16
    plt.subplot(2,2,1)
    plt.plot(times,y, color='blue',LineWidth=plotlinewidth,label="y")
    plt.plot(times,vy, color='green',LineWidth=plotlinewidth,label="vy")
    plt.plot(times,ay, color='black', LineWidth=plotlinewidth,label="ay")
    plt.plot(times,phiseries, color='yellow',LineWidth=plotlinewidth,label="phi")
    plt.legend(loc="best")

    plt.subplot(2,2,2)
    plt.plot(times,z, color='blue',LineWidth=plotlinewidth,label="z")
    plt.plot(times,vz, color='green',LineWidth=plotlinewidth,label="vz")
    plt.plot(times,az, color='black', LineWidth=plotlinewidth,label="az")
    plt.plot(times,thurst, color='yellow',LineWidth=plotlinewidth,label="thurst")
    plt.legend(loc="best")

    plt.subplot(2,2,3)
    plt.plot(-y,z, color='blue',LineWidth=plotlinewidth,label="y-z")
    plt.legend(loc="best")

    plt.subplot(2,2,4)
    plt.plot(times,F1, color='blue',LineWidth=plotlinewidth,label="F1")
    plt.plot(times,F2, color='black', LineWidth=plotlinewidth,label="F2")
    plt.legend(loc="best")

    # print(res)
    # core calculate code
    plt.show()


if __name__ == '__main__':  # 主函数
    main()
