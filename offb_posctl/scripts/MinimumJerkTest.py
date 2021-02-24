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

## core part 0 <<<<
# ay0 vy0 y0 az0 vz0 z0 aytf vytf ytf aztf vztf ztf meshpoint thrustmax angleaccdmax lbz lbv ubv
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
meshpoint=np.linspace(1,0.01,20)
thrustmax=2*9.8
angleaccdmax=25
lbz=0.2
lbv=-5
ubv=5


# Objective
@jit(float64(float64[:]), nopython=True)
def J(x):
    return x[-1]

@jit(float64[:](float64[:]), nopython=True)
def fast_jac(x):
    jac = np.zeros_like(x)
    jac[-1]=1
    # print("jac:",jac)
    return jac


# Constraint
# @jit(float64[:](float64[:]), nopython=True)
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
    # print("1111:",np.vstack((ceq1,ceq2,ceq3,ceq4,ceq5,ceq6)))
    # print("0000:", np.vstack((ceq1,ceq2,ceq3,ceq4,ceq5,ceq6)).shape)
    # print("aztf--------", aztf)
    # print("2222:",np.hstack((ceq1,ceq2,ceq3,ceq4,ceq5,ceq6)))
    # print("3333:",np.vstack((ceq1,ceq2,ceq3,ceq4,ceq5,ceq6)).ravel())
    # return np.hstack((ceq1,ceq2,ceq3,ceq4,ceq5,ceq6))
    #ravel() is used to decrease the dimenions of the array, so the float64[:] may be able to be use. However, it is useless
    return np.hstack((ceq1,ceq2,ceq3,ceq4,ceq5,ceq6)).ravel()


# Constraint
# @jit(float64[:](float64[:]), nopython=True)
def ineqmycon(x):
    global ay0, vy0, y0, az0, vz0, z0, aytf, vytf, ytf, aztf, vztf, ztf, meshpoint, thrustmax, angleaccdmax, lbz, lbv, ubv
    # alpha_y=x[0]
    # beta_y=x[1]
    # gamma_y=x[2]
    #
    # alpha_z=x[3]
    # beta_z=x[4]
    # gamma_z=x[5]
    t=x[0]
    # print("t-----",t)
    tarray=np.array([[60/t**3,-360/t**4,720/t**5],[-24/t**2,168/t**3,-360/t**4],[3/t,-24/t**2,60/t**3]])

    alpha_y,beta_y,gamma_y=np.dot(tarray,np.array([aytf-ay0,vytf-vy0-ay0*t,ytf-y0-vy0*t-0.5*ay0*t**2]))
    alpha_z,beta_z,gamma_z=np.dot(tarray,np.array([aztf-az0,vztf-vz0-az0*t,ztf-z0-vz0*t-0.5*az0*t**2]))
    # print(alpha_z)


    tmesh=t*(np.array(meshpoint))
    angleacc=np.zeros_like(tmesh)
    for i in range(len(tmesh)):
        # print("i===",tmesh)
        t=tmesh[i]
        angleacc[i]=((((alpha_y*t**2)/2 + beta_y*t + gamma_y)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - (((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2)*((2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3))/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)**2 - ((beta_y + alpha_y*t)/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5) - ((beta_z + alpha_z*t)*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 - (2*((alpha_y*t**2)/2 + beta_y*t + gamma_y)*((alpha_z*t**2)/2 + beta_z*t + gamma_z))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + (2*((alpha_z*t**2)/2 + beta_z*t + gamma_z)**2*((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0))/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**3)/(((alpha_y*t**3)/6 + (beta_y*t**2)/2 + gamma_y*t + ay0)**2/((alpha_z*t**3)/6 + (beta_z*t**2)/2 + gamma_z*t + az0 + 49/5)**2 + 1)

    t=x[0]
    thrust=np.sqrt(((alpha_y*tmesh**3)/6 + (beta_y*tmesh**2)/2 + gamma_y*tmesh + ay0)**2 + ((alpha_z*tmesh**3)/6 + (beta_z*tmesh**2)/2 + gamma_z*tmesh + az0 + 49/5)**2)
    c0=t
    # thrust constraints
    c1=2*9.8-thrust
    # print("c1----",c1.shape)
    # z's lower bound  constraints
    c2=-lbz+(alpha_z/120*tmesh**5+beta_z/24*tmesh**4+gamma_z/6*tmesh**3+az0/2*tmesh**2+vz0*tmesh+z0)

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
    # print("--------", t,np.hstack((c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13)))
    # print("--------", (np.hstack((c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13))>0).all())
    return (np.hstack((c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13))>-0.1).all()

def main():

    constraint = [dict(type='ineq', fun=ineqmycon)]

    Initial_guess=np.array([0.1]).ravel()
    lb=-1000
    ub=1000
    mybounds=[(0.5,10)]
    while True:
        Initial_guess=Initial_guess*0.5
        start = time.time()
        # result = minimize(J, Initial_guess, method='SLSQP', jac=fast_jac,tol=1e-4, bounds=mybounds,constraints=constraint)
        while (ineqmycon(Initial_guess)==False):
            Initial_guess=Initial_guess+0.2
            # print(ineqmycon(Initial_guess))
            end = time.time()
            if((end-start)>0.5):
                Initial_guess=0.1
        running_time = end - start
        print('time cost : %.5f sec' % running_time)
        print(Initial_guess[0])
    # times=np.linspace(0,1,100)*Initial_guess
    #
    # t=Initial_guess[0]
    # tarray=np.array([[60/t**3,-360/t**4,720/t**5],[-24/t**2,168/t**3,-360/t**4],[3/t,-24/t**2,60/t**3]])
    #
    # alpha_y,beta_y,gamma_y=np.dot(tarray,np.array([aytf-ay0,vytf-vy0-ay0*t,ytf-y0-vy0*t-0.5*ay0*t**2]))
    # alpha_z,beta_z,gamma_z=np.dot(tarray,np.array([aztf-az0,vztf-vz0-az0*t,ztf-z0-vz0*t-0.5*az0*t**2]))
    # # alpha_y=result.x[0]
    # # beta_y=result.x[1]
    # # gamma_y=result.x[2]
    # #
    # # alpha_z=result.x[3]
    # # beta_z=result.x[4]
    # # gamma_z=result.x[5]
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
