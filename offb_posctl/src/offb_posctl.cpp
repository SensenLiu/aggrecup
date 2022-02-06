/*
 * position_control.cpp
 *
 * Author:mz
 *
 * Time: 2018.11.27
 *
 * 说明: mavros位置控制示例程序
 *      输入：mavros发布的位置/速度信息
 *      输出：无人机的推力和姿态信息
 *      采用位置环/速度环串级PID控制，位置环P控制，速度环PID控制
 */
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <queue>
#include <vector>
#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <Eigen/Eigen>
#include <Eigen/Geometry> 
#include <Eigen/Core> 


#include <ros/ros.h>
#include "Parameter.h"
#include <PID.h>
#include <FILTER.h>
#include <DOB.h>


//topic
#include <geometry_msgs/Point.h>
#include <geometry_msgs/Vector3.h>
#include <geometry_msgs/Quaternion.h>
#include <geometry_msgs/Pose.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <mavros_msgs/State.h>
#include <mavros_msgs/AttitudeTarget.h>
#include <mavros_msgs/CommandBool.h>
#include <mavros_msgs/SetMode.h>
#include <mavros_msgs/PositionTarget.h>
#include <nav_msgs/Odometry.h>
#include <std_msgs/Bool.h>
#include <std_msgs/Float32.h>
#include <std_msgs/Float64.h>
#include <mavros_msgs/Thrust.h>
#include <mavros_msgs/AttitudeTarget.h>
#include <sensor_msgs/Imu.h>
#include <nav_msgs/Path.h>

#include "ros/ros.h"
#include "std_msgs/Float32.h"
#include <chrono>
#include "iomanip"
#include <thread>
#include <unistd.h>
#include <mutex>
#include "offb_posctl/controlstate.h"
#include "nlopt.hpp"
using namespace Eigen;//释放eigen命名空间 矩阵库
using namespace std;

///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>全 局 变 量<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

mavros_msgs::State current_state;           //无人机当前状态
nav_msgs::Odometry pose_drone_odom;       //读入的无人机drone当前位置，x，y，z+姿态
sensor_msgs::Imu pose_drone_Imu;       //读入的无人机drone当前位置，x，y，z+姿态
nav_msgs::Odometry pose_car_odom;       //读入car当前位置
nav_msgs::Odometry wallpostwist_odom;       //读入car当前位置
geometry_msgs::TwistStamped vel_drone;      //读入的无人机当前速度 线速度+角速度
geometry_msgs::Quaternion orientation_target;   //发给无人机的姿态指令  四元数
geometry_msgs::Vector3 angle_target;   //欧拉角
geometry_msgs::Vector3 vel_read;   //期望速度
geometry_msgs::Point plane_expected_position; //车的零点和飞机的零点差3m，根据车的当前位置计算飞机位置
geometry_msgs::Point plane_expected_velocity; //车的零点和飞机的零点差3m，根据车的当前位置计算飞机位置
geometry_msgs::Point plane_expected_acceleration; //车的零点和飞机的零点差3m，根据车的当前位置计算飞机位置

std_msgs::Float64 plane_real_alt; //control前
mavros_msgs::AttitudeTarget target_atti_thrust_msg; //最终发布的消息 油门+角度
nav_msgs::Odometry  planned_postwist_msg;

nav_msgs::Odometry current_relativepostwist_msg;
geometry_msgs::Vector3 targeterror_msg;
geometry_msgs::Vector3 temp_angle;
geometry_msgs::Vector3 rpy;
offb_posctl::controlstate controlstatearray_msg;
geometry_msgs::Vector3 angle_receive;       //读入的无人机姿态（欧拉角）
geometry_msgs::PoseStamped  pos_drone;      //读入的无人机当前位置
geometry_msgs::PoseStamped  pos_wall; //读入的无人机上一次位置
geometry_msgs::TwistStamped vel_wall;      //读入的无人机当前速度
float thrust_target;        //期望推力
float Yaw_Init;
float Yaw_Locked = 0;           //锁定的偏航角(一般锁定为0)
bool got_initial_point = false;
PID PIDVX, PIDVY, PIDVZ;    //声明PID类
Parameter param;
std::ofstream logfile;

///for psopt
float px_ini = -0.0;
float pz_ini = 0;
float py_ini=0;
float vx_ini = -0.0;
float vz_ini = 0.0;
float vy_ini=0;
float ax_ini=0;
float ay_ini=0;
float az_ini = 0.0;

float ax_des=0;
float ay_des=0;
float az_des=0.0;

float phi_ini=0;
float thrust_ini=9.8;

int controlfreq=30;
int discretizedpointpersecond = controlfreq;
int controlcounter=0;
int quad_state=0;//0 climbhover, 1 AggressiveFly, 2 AdhesionPhase, 3 keep current state hover, 4 AdhesionSuccess
FILTER vy_Der(20),vz_Der(20),ay_Filter(20),az_Filter(20);
FILTER ax_sgFilter(51,11,2),ay_sgFilter(51,11,2),az_sgFilter(51,11,2);
DOB x_dob(5,0.2),y_dob(5,0.2),z_dob(5,0.2);
double x_acccomp=0,y_acccomp=0,z_acccomp=0;
///for pitch compensation
float comp_integrate,comp_last_error;
float comp_kp = 0,comp_ki = 0,comp_kd = 0;

bool planstopflag=false;
bool startattitudecotrolflag=false;
bool anglereachedflag=false;

float feedforwardcoefx=1.0,feedforwardcoefy=1.0,feedforwardcoefz=1.0;
float offsetdrone_wall_x=10000;
float offsetdrone_wall_y=1.72-1.87;
float offsetdrone_wall_z=1.39-1.52;//phi 1.22


///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>terminal replanning<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
int timeswitchnumber = 4;
typedef struct {
    double T, k, a0, u1, u2, x10, x20, x30, x40, x1tf, x2tf, x3tf, x4tf, lmd1, lmd2, lmd3, lmd4;
} my_cost_data;
my_cost_data costdata = {0.6, 1.5, 0.3, -1, 1, 8.7, 4.2, 0.96, 1.03, 11, 4.2, 1.02, 1.03, 5, 1, 5, 1};
nlopt::opt opt(nlopt::LD_SLSQP, timeswitchnumber);
bool targetangleachievedflag= false;
int ocp_startindex=0;
Eigen::ArrayXf ocp_thrust_array=ArrayXf::Ones(1)*9.8*max(0.3,cos(1.22));;
double ytf_ocp=0,vytf_ocp=0,ztf_ocp=0,vztf_ocp=0;
///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>声 明 函 数<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

//欧拉角转四元数
geometry_msgs::Quaternion euler2quaternion(float roll, float pitch, float yaw);//geometry_msgs的Quaternion类型的函数
geometry_msgs::Vector3 quaternion2euler(float x, float y, float z, float w);

float get_ros_time(ros::Time time_begin);                                            //获取ros当前时间
int pix_controller(float cur_time);
void vector3dLimit(Vector3d &v, double limit) ; ///limit should be positive
Vector3d vectorElementMultiply(Vector3d v1, Vector3d v2);
void data_log(std::ofstream &logfile, float cur_time);
float pitch_compensation(float theta_dc, float theta_dn, float yk[],float uk[]);

///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>回 调 函 数<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

void state_cb(const mavros_msgs::State::ConstPtr &msg){
    current_state = *msg;

}//当有消息到达topic时会自动调用一次
bool planeupdateflag= false;
void plane_pos_cb(const nav_msgs::Odometry::ConstPtr &msg){
    pose_drone_odom = *msg;//pose_drone_odom is nav_msgs::Odometry type
    planeupdateflag= true;

    pos_drone.header=pose_drone_odom.header;
    pos_drone.pose=pose_drone_odom.pose.pose;
}
void plane_imu_cb(const sensor_msgs::Imu::ConstPtr &msg){

    pose_drone_Imu = *msg;
    angle_receive = quaternion2euler(pose_drone_Imu.orientation.x, pose_drone_Imu.orientation.y, pose_drone_Imu.orientation.z, pose_drone_Imu.orientation.w);
//    ROS_ERROR_STREAM_THROTTLE(1,"drone angle_receive.x:"<<angle_receive.x);
    Eigen::Quaterniond yawRot_q=Eigen::AngleAxisd(-angle_receive.z, Vector3d::UnitZ()) * Eigen::AngleAxisd(0, Vector3d::UnitY()) * Eigen::AngleAxisd(0, Vector3d::UnitX());
    Eigen::Quaterniond current_q(pose_drone_Imu.orientation.w, pose_drone_Imu.orientation.x, pose_drone_Imu.orientation.y, pose_drone_Imu.orientation.z);
    Eigen::Quaterniond current_acc(0, pose_drone_Imu.linear_acceleration.x, pose_drone_Imu.linear_acceleration.y, pose_drone_Imu.linear_acceleration.z);
    Eigen::Quaterniond accinword=yawRot_q*current_q*current_acc*current_q.conjugate()*yawRot_q.conjugate();

    ax_ini=accinword.x();
    ay_ini=accinword.y();
    az_ini=accinword.z()-9.8;

    Eigen::Quaterniond des_acc(0, 0, 0, target_atti_thrust_msg.thrust/param.THR_HOVER*9.8);
    Eigen::Quaterniond desaccinword=yawRot_q*current_q*des_acc*current_q.conjugate()*yawRot_q.conjugate();
    ax_des=desaccinword.x();
    ay_des=desaccinword.y();
    az_des=desaccinword.z()-9.8;
}

void plane_vel_cb(const geometry_msgs::TwistStamped::ConstPtr &msg){
    vel_drone = *msg;

    vel_read=vel_drone.twist.linear;
}

void car_pos_cb(const nav_msgs::Odometry::ConstPtr &msg) {
    pose_car_odom = *msg;

//    plane_expected_position.x = pose_car_odom.pose.pose.position.x; //-1 means axis difference
//    plane_expected_position.y = pose_car_odom.pose.pose.position.y; //-1 means axis difference
//    plane_expected_position.z = pose_car_odom.pose.pose.position.z + 0.5;
}
bool wallstate_receivedflag= false;
void wall_twist_cb(const nav_msgs::Odometry::ConstPtr &msg) {
    wallpostwist_odom = *msg;
    wallstate_receivedflag= true;
//    plane_expected_position.x = pose_car_odom.pose.pose.position.x; //-1 means axis difference
//    plane_expected_position.y = pose_car_odom.pose.pose.position.y; //-1 means axis difference
//    plane_expected_position.z = pose_car_odom.pose.pose.position.z + 0.5;

    pos_wall.header=wallpostwist_odom.header;
    pos_wall.pose=wallpostwist_odom.pose.pose;
    vel_wall.header=wallpostwist_odom.header;
    vel_wall.twist=wallpostwist_odom.twist.twist;
}
void plane_alt_cb(const std_msgs::Float64::ConstPtr &msg){
    plane_real_alt = *msg;
}
bool contstaterecieveflag= false;
void controlstate_cb(const offb_posctl::controlstate::ConstPtr &msg)
{
    if(planstopflag== false)
    {
        controlstatearray_msg = *msg;
        controlcounter = controlstatearray_msg.inicounter;
        if(contstaterecieveflag == false)//第一次回调时初始化，之后这个flag一直是true
        {
            contstaterecieveflag= true;
            discretizedpointpersecond=controlstatearray_msg.discrepointpersecond;
        }
    }

}

///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>last segment thrust ocp<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    my_cost_data *d = reinterpret_cast<my_cost_data*>(my_func_data);
    double T=d->T, k=d->k, a0=d->a0, u1=d->u1, u2=d->u2, x10=d->x10, x20=d->x20, x30=d->x30, x40=d->x40, x1tf=d->x1tf,
            x2tf=d->x2tf, x3tf=d->x3tf, x4tf=d->x4tf, lmd1=d->lmd1, lmd2=d->lmd2, lmd3=d->lmd3, lmd4=d->lmd4;
    double t1=x[0],t2=x[1],t3=x[2],t4=x[3];
    if (!grad.empty()) {
        grad[0] = 2*lmd2*(u1*sin(a0 + k*t1) - u2*sin(a0 + k*t1))*(x2tf - x20 - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/k + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/k + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/k + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/k - (u1*(cos(a0 + k*t1) - cos(a0)))/k) - 2*lmd3*((k*u1*sin(a0) - k*u1*sin(a0 + k*t1))/pow(k,2) - u2*cos(a0 + k*t1)*(t1 - t2))*(x30 - x3tf + T*x40 - (u1*(cos(a0 + k*t1) - cos(a0)) + k*t1*u1*sin(a0))/pow(k,2) - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/pow(k,2) + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/pow(k,2) + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/pow(k,2) + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/pow(k,2) - (u1*sin(a0 + k*t4)*(T - t4))/k + (u2*sin(a0 + k*t1)*(t1 - t2))/k + (u1*sin(a0 + k*t2)*(t2 - t3))/k + (u2*sin(a0 + k*t3)*(t3 - t4))/k) - 2*lmd1*((k*u1*cos(a0) - k*u1*cos(a0 + k*t1))/pow(k,2) + u2*sin(a0 + k*t1)*(t1 - t2))*(x10 - x1tf + T*x20 + (u1*(sin(a0 + k*t1) - sin(a0)) - k*t1*u1*cos(a0))/pow(k,2) + (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/pow(k,2) - (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/pow(k,2) - (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/pow(k,2) - (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/pow(k,2) - (u1*cos(a0 + k*t4)*(T - t4))/k + (u2*cos(a0 + k*t1)*(t1 - t2))/k + (u1*cos(a0 + k*t2)*(t2 - t3))/k + (u2*cos(a0 + k*t3)*(t3 - t4))/k) - 2*lmd4*(u1*cos(a0 + k*t1) - u2*cos(a0 + k*t1))*(x4tf - x40 - (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/k + (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/k + (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/k + (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/k - (u1*(sin(a0 + k*t1) - sin(a0)))/k);
        grad[1] = 2*lmd4*(u1*cos(a0 + k*t2) - u2*cos(a0 + k*t2))*(x4tf - x40 - (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/k + (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/k + (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/k + (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/k - (u1*(sin(a0 + k*t1) - sin(a0)))/k) - 2*lmd2*(u1*sin(a0 + k*t2) - u2*sin(a0 + k*t2))*(x2tf - x20 - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/k + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/k + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/k + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/k - (u1*(cos(a0 + k*t1) - cos(a0)))/k) - 2*lmd1*((u2*cos(a0 + k*t1))/k - (u2*cos(a0 + k*t2))/k + u1*sin(a0 + k*t2)*(t2 - t3))*(x10 - x1tf + T*x20 + (u1*(sin(a0 + k*t1) - sin(a0)) - k*t1*u1*cos(a0))/pow(k,2) + (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/pow(k,2) - (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/pow(k,2) - (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/pow(k,2) - (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/pow(k,2) - (u1*cos(a0 + k*t4)*(T - t4))/k + (u2*cos(a0 + k*t1)*(t1 - t2))/k + (u1*cos(a0 + k*t2)*(t2 - t3))/k + (u2*cos(a0 + k*t3)*(t3 - t4))/k) + 2*lmd3*((u2*sin(a0 + k*t2))/k - (u2*sin(a0 + k*t1))/k + u1*cos(a0 + k*t2)*(t2 - t3))*(x30 - x3tf + T*x40 - (u1*(cos(a0 + k*t1) - cos(a0)) + k*t1*u1*sin(a0))/pow(k,2) - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/pow(k,2) + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/pow(k,2) + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/pow(k,2) + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/pow(k,2) - (u1*sin(a0 + k*t4)*(T - t4))/k + (u2*sin(a0 + k*t1)*(t1 - t2))/k + (u1*sin(a0 + k*t2)*(t2 - t3))/k + (u2*sin(a0 + k*t3)*(t3 - t4))/k);
        grad[2] = 2*lmd2*(u1*sin(a0 + k*t3) - u2*sin(a0 + k*t3))*(x2tf - x20 - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/k + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/k + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/k + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/k - (u1*(cos(a0 + k*t1) - cos(a0)))/k) - 2*lmd1*((u1*cos(a0 + k*t2))/k - (u1*cos(a0 + k*t3))/k + u2*sin(a0 + k*t3)*(t3 - t4))*(x10 - x1tf + T*x20 + (u1*(sin(a0 + k*t1) - sin(a0)) - k*t1*u1*cos(a0))/pow(k,2) + (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/pow(k,2) - (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/pow(k,2) - (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/pow(k,2) - (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/pow(k,2) - (u1*cos(a0 + k*t4)*(T - t4))/k + (u2*cos(a0 + k*t1)*(t1 - t2))/k + (u1*cos(a0 + k*t2)*(t2 - t3))/k + (u2*cos(a0 + k*t3)*(t3 - t4))/k) - 2*lmd4*(u1*cos(a0 + k*t3) - u2*cos(a0 + k*t3))*(x4tf - x40 - (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/k + (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/k + (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/k + (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/k - (u1*(sin(a0 + k*t1) - sin(a0)))/k) + 2*lmd3*((u1*sin(a0 + k*t3))/k - (u1*sin(a0 + k*t2))/k + u2*cos(a0 + k*t3)*(t3 - t4))*(x30 - x3tf + T*x40 - (u1*(cos(a0 + k*t1) - cos(a0)) + k*t1*u1*sin(a0))/pow(k,2) - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/pow(k,2) + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/pow(k,2) + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/pow(k,2) + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/pow(k,2) - (u1*sin(a0 + k*t4)*(T - t4))/k + (u2*sin(a0 + k*t1)*(t1 - t2))/k + (u1*sin(a0 + k*t2)*(t2 - t3))/k + (u2*sin(a0 + k*t3)*(t3 - t4))/k);
        grad[3] = 2*lmd1*((u2*cos(a0 + k*t4))/k - (u2*cos(a0 + k*t3))/k + u1*sin(a0 + k*t4)*(T - t4))*(x10 - x1tf + T*x20 + (u1*(sin(a0 + k*t1) - sin(a0)) - k*t1*u1*cos(a0))/pow(k,2) + (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/pow(k,2) - (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/pow(k,2) - (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/pow(k,2) - (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/pow(k,2) - (u1*cos(a0 + k*t4)*(T - t4))/k + (u2*cos(a0 + k*t1)*(t1 - t2))/k + (u1*cos(a0 + k*t2)*(t2 - t3))/k + (u2*cos(a0 + k*t3)*(t3 - t4))/k) - 2*lmd2*(u1*sin(a0 + k*t4) - u2*sin(a0 + k*t4))*(x2tf - x20 - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/k + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/k + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/k + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/k - (u1*(cos(a0 + k*t1) - cos(a0)))/k) + 2*lmd4*(u1*cos(a0 + k*t4) - u2*cos(a0 + k*t4))*(x4tf - x40 - (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/k + (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/k + (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/k + (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/k - (u1*(sin(a0 + k*t1) - sin(a0)))/k) - 2*lmd3*((u2*sin(a0 + k*t3))/k - (u2*sin(a0 + k*t4))/k + u1*cos(a0 + k*t4)*(T - t4))*(x30 - x3tf + T*x40 - (u1*(cos(a0 + k*t1) - cos(a0)) + k*t1*u1*sin(a0))/pow(k,2) - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/pow(k,2) + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/pow(k,2) + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/pow(k,2) + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/pow(k,2) - (u1*sin(a0 + k*t4)*(T - t4))/k + (u2*sin(a0 + k*t1)*(t1 - t2))/k + (u1*sin(a0 + k*t2)*(t2 - t3))/k + (u2*sin(a0 + k*t3)*(t3 - t4))/k);
    }
//    double J= lmd1*pow((x10 - x1tf + T*x20 + (u1*(sin(a0 + k*t1) - sin(a0)) - k*t1*u1*cos(a0))/pow(k,2) + (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/pow(k,2) - (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/pow(k,2) - (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/pow(k,2) - (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/pow(k,2) - (u1*cos(a0 + k*t4)*(T - t4))/k + (u2*cos(a0 + k*t1)*(t1 - t2))/k + (u1*cos(a0 + k*t2)*(t2 - t3))/k + (u2*cos(a0 + k*t3)*(t3 - t4))/k),2) + lmd4*pow((x4tf - x40 - (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/k + (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/k + (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/k + (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/k - (u1*(sin(a0 + k*t1) - sin(a0)))/k),2) + lmd2*pow((x2tf - x20 - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/k + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/k + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/k + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/k - (u1*(cos(a0 + k*t1) - cos(a0)))/k),2) + lmd3*pow((x30 - x3tf + T*x40 - (u1*(cos(a0 + k*t1) - cos(a0)) + k*t1*u1*sin(a0))/pow(k,2) - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/pow(k,2) + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/pow(k,2) + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/pow(k,2) + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/pow(k,2) - (u1*sin(a0 + k*t4)*(T - t4))/k + (u2*sin(a0 + k*t1)*(t1 - t2))/k + (u1*sin(a0 + k*t2)*(t2 - t3))/k + (u2*sin(a0 + k*t3)*(t3 - t4))/k),2);
//
//    cout<<"cost J: "<<J<<endl;

    return  lmd1*pow((x10 - x1tf + T*x20 + (u1*(sin(a0 + k*t1) - sin(a0)) - k*t1*u1*cos(a0))/pow(k,2) + (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/pow(k,2) - (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/pow(k,2) - (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/pow(k,2) - (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/pow(k,2) - (u1*cos(a0 + k*t4)*(T - t4))/k + (u2*cos(a0 + k*t1)*(t1 - t2))/k + (u1*cos(a0 + k*t2)*(t2 - t3))/k + (u2*cos(a0 + k*t3)*(t3 - t4))/k),2) + lmd4*pow((x4tf - x40 - (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/k + (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/k + (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/k + (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/k - (u1*(sin(a0 + k*t1) - sin(a0)))/k),2) + lmd2*pow((x2tf - x20 - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/k + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/k + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/k + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/k - (u1*(cos(a0 + k*t1) - cos(a0)))/k),2) + lmd3*pow((x30 - x3tf + T*x40 - (u1*(cos(a0 + k*t1) - cos(a0)) + k*t1*u1*sin(a0))/pow(k,2) - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/pow(k,2) + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/pow(k,2) + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/pow(k,2) + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/pow(k,2) - (u1*sin(a0 + k*t4)*(T - t4))/k + (u2*sin(a0 + k*t1)*(t1 - t2))/k + (u1*sin(a0 + k*t2)*(t2 - t3))/k + (u2*sin(a0 + k*t3)*(t3 - t4))/k),2);
}


void multiconstraint(unsigned m, double *result, unsigned n, const double *x, double *gradient, void *func_data) // ref:https://blog.csdn.net/weixin_33937499/article/details/93488349 and https://nlopt.readthedocs.io/en/latest/NLopt_Reference/
{
    if(gradient){
        gradient[0]= 1.0;
        gradient[1]= -1.0;
        gradient[2]= 0.0;
        gradient[3]= 0.0;

        gradient[4]= 0.0;
        gradient[5]= 1.0;
        gradient[6]= -1.0;
        gradient[7]= 0.0;

        gradient[8]= 0.0;
        gradient[9]= 0.0;
        gradient[10]= 1.0;
        gradient[11]= -1.0;
    }
    result[0] = x[0]-x[1];
    result[1] = x[1]-x[2];
    result[2] = x[2]-x[3];
}

bool terminalreplanning()
{
    if(targetangleachievedflag==false && controlstatearray_msg.arraylength)
    {
//        costdata = {0.6, 1.5, 0.3, -1, 1, 8.7, 4.2, 0.96, 1.03, 11, 4.2, 1.02, 1.03, 1, 1, 1, 1};
        double currenttime=ros::Time::now().toSec();
        int startindex=round((currenttime-controlstatearray_msg.header.stamp.toSec())*controlstatearray_msg.discrepointpersecond)+controlstatearray_msg.inicounter <=controlstatearray_msg.arraylength ? round((currenttime-controlstatearray_msg.header.stamp.toSec())*controlstatearray_msg.discrepointpersecond)+controlstatearray_msg.inicounter : controlstatearray_msg.arraylength;
        double plannedend_roll=min(-atan2(controlstatearray_msg.stateAYarray[controlstatearray_msg.arraylength-1],(controlstatearray_msg.stateAZarray[controlstatearray_msg.arraylength-1]+9.8)),1.57);
        int array_length=0;
        double angleoffset=0.2;
        costdata.a0=angle_receive.x;
        if(controlstatearray_msg.arraylength-startindex>3)
        {
            double plannnedstart_roll=-atan2(controlstatearray_msg.stateAYarray[startindex],(controlstatearray_msg.stateAZarray[startindex]+9.8));
            costdata.k=(plannedend_roll-plannnedstart_roll)/((double)controlstatearray_msg.arraylength/(double)controlstatearray_msg.discrepointpersecond);
        }

        if (costdata.k>=0)
        {
            costdata.T=(plannedend_roll-angleoffset-costdata.a0)/costdata.k;
            if(costdata.T<0.03)
            {
                targetangleachievedflag= true;
                ocp_thrust_array=ArrayXf::Ones(1)*9.8*max(0.3,cos(plannedend_roll));
                ocp_startindex=0;
//                ROS_WARN_STREAM( "---------------------targetangleachievedflag---------------"<<targetangleachievedflag<<" 9.8*max(0.3,cos(plannedend_roll):"<<9.8*max(0.3,cos(plannedend_roll)));
//                cout<<" ocp_thrust_array:"<<ocp_thrust_array<<endl;
                return false;
            }
            array_length=round(costdata.T*controlstatearray_msg.discrepointpersecond);
        }

        Eigen::ArrayXf thrust_array = ArrayXf::Ones(array_length)*9.8*max(0.3,cos(plannedend_roll));
        Eigen::ArrayXf compensatethrust=ArrayXf::Zero(array_length);
        Eigen::ArrayXf roll_array = ArrayXf::LinSpaced(array_length,costdata.a0,plannedend_roll);
        Eigen::ArrayXf thrustsin_array = ArrayXf::Zero(array_length);
        Eigen::ArrayXf thrustcos_array = ArrayXf::Zero(array_length);
        Eigen::MatrixXf interalmatrix=MatrixXf ::Ones(array_length,array_length).triangularView<Eigen::Upper>();
        for(int i=startindex;i<controlstatearray_msg.arraylength;i++)
        {
            if((i-startindex)>=array_length)
            {
                break;
            }
            thrust_array(i-startindex)=sqrt(pow(controlstatearray_msg.stateAYarray[i],2)+pow(controlstatearray_msg.stateAZarray[i]+9.8,2));
        }

        ocp_thrust_array=thrust_array+compensatethrust;
        thrustsin_array=thrust_array*sin(roll_array);
        thrustsin_array=(thrustsin_array.matrix()*interalmatrix*1.0/controlstatearray_msg.discrepointpersecond).array();
        thrustcos_array=thrust_array*cos(roll_array);
        thrustcos_array=(thrustcos_array.matrix()*interalmatrix*1.0/controlstatearray_msg.discrepointpersecond).array();

        double lowthrust=min(0.3*9.8-thrust_array.minCoeff(),0.0);
        double upperthrust=max(2*9.8-thrust_array.maxCoeff(),0.0);
//        cout<<"thrust_array.minCoeff(): "<<thrust_array.minCoeff()<<" thrust_array.maxCoeff():"<<thrust_array.maxCoeff()<<endl;
        ytf_ocp=(currenttime+costdata.T-(controlstatearray_msg.header.stamp.toSec()+(double)controlstatearray_msg.arraylength/controlstatearray_msg.discrepointpersecond))*controlstatearray_msg.rendezvouswall_vy+controlstatearray_msg.rendezvouswall_y;
        vytf_ocp=controlstatearray_msg.stateVYarray[controlstatearray_msg.arraylength-1];
        ztf_ocp=(currenttime+costdata.T-(controlstatearray_msg.header.stamp.toSec()+(double)controlstatearray_msg.arraylength/controlstatearray_msg.discrepointpersecond))*controlstatearray_msg.rendezvouswall_vz+controlstatearray_msg.rendezvouswall_z;
        vztf_ocp=controlstatearray_msg.stateVZarray[controlstatearray_msg.arraylength-1];

        costdata.u1=lowthrust;costdata.u2=upperthrust;
        costdata.x10=current_relativepostwist_msg.pose.pose.position.y;
        costdata.x20=current_relativepostwist_msg.twist.twist.linear.y;
        costdata.x30=current_relativepostwist_msg.pose.pose.position.z;
        costdata.x40=current_relativepostwist_msg.twist.twist.linear.z;

        costdata.x1tf=ytf_ocp+thrustsin_array.sum()*1.0/controlstatearray_msg.discrepointpersecond;
        costdata.x2tf=vytf_ocp+thrustsin_array(array_length-1);
        costdata.x3tf=ztf_ocp-thrustcos_array.sum()*1.0/controlstatearray_msg.discrepointpersecond+9.8*costdata.T*costdata.T/2;
        costdata.x4tf=vztf_ocp-thrustcos_array(array_length-1)+9.8*costdata.T;

        std::vector<double> ub(timeswitchnumber, costdata.T);
        opt.set_upper_bounds(ub);

        std::vector<double> x1(timeswitchnumber);
        x1[0] = 0.0 / (timeswitchnumber + 1) * costdata.T;
        x1[1] = 2.0 / (timeswitchnumber + 1) * costdata.T;
        x1[2] = 3.0 / (timeswitchnumber + 1) * costdata.T;
        x1[3] = 4.0 / (timeswitchnumber + 1) * costdata.T;
        std::vector<double> x2{x1[0],x1[1],x1[2],x1[3]};
        double minf1=+HUGE_VAL,minf2=+HUGE_VAL;
        opt.set_min_objective(myvfunc, &costdata);
        try {
            nlopt::result result = opt.optimize(x1, minf1);
            std::cout << "found minimum at f(" << x1[0] << "," << x1[1] << "," << x1[2] << "," << x1[3] << ") = "
                      << std::setprecision(3) << minf1 << "remaining T "<<costdata.T<<std::endl;
        }
        catch (std::exception &e) {
            std::cout << "nlopt failed: " << e.what() << std::endl;
        }
        costdata.u1=upperthrust;costdata.u2=lowthrust;
        opt.set_min_objective(myvfunc, &costdata);
        try {
            nlopt::result result = opt.optimize(x2, minf2);
            std::cout << "found minimum at f(" << x2[0] << "," << x2[1] << "," << x2[2] << "," << x2[3] << ") = "
                      << std::setprecision(3) << minf2 << "remaining T "<<costdata.T<<std::endl;
        }
        catch (std::exception &e) {
            std::cout << "nlopt failed: " << e.what() << std::endl;
        }
        if(minf1==+HUGE_VAL&&minf2==+HUGE_VAL)
        {
            return false;
        }
        if(minf1<minf2) {
            Eigen::ArrayXf tempcompensatethrust=ArrayXf::Zero(array_length+1);// +1 is to avoid the index over in switchMatrix.block(1,array_length,p,q) block operation
            Eigen::Array4i timenodevector(round(x1[0]*controlstatearray_msg.discrepointpersecond),round(x1[1]*controlstatearray_msg.discrepointpersecond),round(x1[2]*controlstatearray_msg.discrepointpersecond),round(x1[3]*controlstatearray_msg.discrepointpersecond));
            tempcompensatethrust.segment(0,timenodevector(0))=tempcompensatethrust.segment(0,timenodevector(0))+lowthrust;
            tempcompensatethrust.segment(timenodevector(0),timenodevector(1)-timenodevector(0))=tempcompensatethrust.segment(timenodevector(0),timenodevector(1)-timenodevector(0))+upperthrust;
            tempcompensatethrust.segment(timenodevector(1),timenodevector(2)-timenodevector(1))=tempcompensatethrust.segment(timenodevector(1),timenodevector(2)-timenodevector(1))+lowthrust;
            tempcompensatethrust.segment(timenodevector(2),timenodevector(3)-timenodevector(2))=tempcompensatethrust.segment(timenodevector(2),timenodevector(3)-timenodevector(2))+upperthrust;
            tempcompensatethrust.segment(timenodevector(3),array_length-timenodevector(3))=tempcompensatethrust.segment(timenodevector(3),array_length-timenodevector(3))+lowthrust;
            compensatethrust=tempcompensatethrust.segment(0,array_length);

            double T=costdata.T, k=costdata.k, a0=costdata.a0, u1=costdata.u1, u2=costdata.u2, x10=costdata.x10, x20=costdata.x20, x30=costdata.x30,
            x40=costdata.x40, x1tf=costdata.x1tf, x2tf=costdata.x2tf, x3tf=costdata.x3tf, x4tf=costdata.x4tf, lmd1=costdata.lmd1, lmd2=costdata.lmd2, lmd3=costdata.lmd3, lmd4=costdata.lmd4;
            double t1=x1[0],t2=x1[1],t3=x1[2],t4=x1[3];
            double y_error=lmd1*pow((x10 - x1tf + T*x20 + (u1*(sin(a0 + k*t1) - sin(a0)) - k*t1*u1*cos(a0))/pow(k,2) + (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/pow(k,2) - (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/pow(k,2) - (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/pow(k,2) - (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/pow(k,2) - (u1*cos(a0 + k*t4)*(T - t4))/k + (u2*cos(a0 + k*t1)*(t1 - t2))/k + (u1*cos(a0 + k*t2)*(t2 - t3))/k + (u2*cos(a0 + k*t3)*(t3 - t4))/k),2);
            double vy_error=lmd2*pow((x2tf - x20 - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/k + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/k + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/k + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/k - (u1*(cos(a0 + k*t1) - cos(a0)))/k),2);
            double z_error=lmd3*pow((x30 - x3tf + T*x40 - (u1*(cos(a0 + k*t1) - cos(a0)) + k*t1*u1*sin(a0))/pow(k,2) - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/pow(k,2) + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/pow(k,2) + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/pow(k,2) + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/pow(k,2) - (u1*sin(a0 + k*t4)*(T - t4))/k + (u2*sin(a0 + k*t1)*(t1 - t2))/k + (u1*sin(a0 + k*t2)*(t2 - t3))/k + (u2*sin(a0 + k*t3)*(t3 - t4))/k),2);
            double vz_error=lmd4*pow((x4tf - x40 - (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/k + (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/k + (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/k + (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/k - (u1*(sin(a0 + k*t1) - sin(a0)))/k),2);
            if(y_error+z_error)
            {
                costdata.lmd1=y_error/(y_error+z_error);
                costdata.lmd3=10*z_error/(y_error+z_error);
            }
            if(vy_error+vz_error)
            {
                costdata.lmd2=vy_error/(vy_error+vz_error);
                costdata.lmd4=10*vz_error/(vy_error+vz_error);
            }

            //            Eigen::MatrixXf switchMatrix=MatrixXf::Zero(2,array_length+1);// +1 is to avoid the index over in switchMatrix.block(1,array_length,p,q) block operation
//            Eigen::Vector2f inputset(lowthrust,upperthrust);
//            switchMatrix.block(0,0,1,timenodevector(0))=Eigen::MatrixXf::Ones(1,timenodevector(0));
////            if(timenodevector(1)-timenodevector(0))
//            {
//                switchMatrix.block(1,timenodevector(0),1,timenodevector(1)-timenodevector(0))=Eigen::MatrixXf::Ones(1,timenodevector(1)-timenodevector(0));
//            }
//            switchMatrix.block(0,timenodevector(1),1,timenodevector(2)-timenodevector(1))=Eigen::MatrixXf::Ones(1,timenodevector(2)-timenodevector(1));
//            switchMatrix.block(1,timenodevector(2),1,timenodevector(3)-timenodevector(2))=Eigen::MatrixXf::Ones(1,timenodevector(3)-timenodevector(2));
//            switchMatrix.block(0,timenodevector(3),1,array_length-timenodevector(3))=Eigen::MatrixXf::Ones(1,array_length-timenodevector(3));
//            compensatethrust=(inputset * switchMatrix).array().block(0,0,1,array_length);
        } else{
            Eigen::ArrayXf tempcompensatethrust=ArrayXf::Zero(array_length+1);
            Eigen::Array4i timenodevector(round(x2[0]*controlstatearray_msg.discrepointpersecond),round(x2[1]*controlstatearray_msg.discrepointpersecond),round(x2[2]*controlstatearray_msg.discrepointpersecond),round(x2[3]*controlstatearray_msg.discrepointpersecond));

            tempcompensatethrust.segment(0,timenodevector(0))=tempcompensatethrust.segment(0,timenodevector(0))+upperthrust;
            tempcompensatethrust.segment(timenodevector(0),timenodevector(1)-timenodevector(0))=tempcompensatethrust.segment(timenodevector(0),timenodevector(1)-timenodevector(0))+lowthrust;
            tempcompensatethrust.segment(timenodevector(1),timenodevector(2)-timenodevector(1))=tempcompensatethrust.segment(timenodevector(1),timenodevector(2)-timenodevector(1))+upperthrust;
            tempcompensatethrust.segment(timenodevector(2),timenodevector(3)-timenodevector(2))=tempcompensatethrust.segment(timenodevector(2),timenodevector(3)-timenodevector(2))+lowthrust;
            tempcompensatethrust.segment(timenodevector(3),array_length-timenodevector(3))=tempcompensatethrust.segment(timenodevector(3),array_length-timenodevector(3))+upperthrust;
            compensatethrust=tempcompensatethrust.segment(0,array_length);

            double T=costdata.T, k=costdata.k, a0=costdata.a0, u1=costdata.u1, u2=costdata.u2, x10=costdata.x10, x20=costdata.x20, x30=costdata.x30,
                    x40=costdata.x40, x1tf=costdata.x1tf, x2tf=costdata.x2tf, x3tf=costdata.x3tf, x4tf=costdata.x4tf, lmd1=costdata.lmd1, lmd2=costdata.lmd2, lmd3=costdata.lmd3, lmd4=costdata.lmd4;
            double t1=x2[0],t2=x2[1],t3=x2[2],t4=x2[3];
            double y_error=lmd1*pow((x10 - x1tf + T*x20 + (u1*(sin(a0 + k*t1) - sin(a0)) - k*t1*u1*cos(a0))/pow(k,2) + (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/pow(k,2) - (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/pow(k,2) - (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/pow(k,2) - (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/pow(k,2) - (u1*cos(a0 + k*t4)*(T - t4))/k + (u2*cos(a0 + k*t1)*(t1 - t2))/k + (u1*cos(a0 + k*t2)*(t2 - t3))/k + (u2*cos(a0 + k*t3)*(t3 - t4))/k),2);
            double vy_error=lmd2*pow((x2tf - x20 - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/k + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/k + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/k + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/k - (u1*(cos(a0 + k*t1) - cos(a0)))/k),2);
            double z_error=lmd3*pow((x30 - x3tf + T*x40 - (u1*(cos(a0 + k*t1) - cos(a0)) + k*t1*u1*sin(a0))/pow(k,2) - (u1*(cos(a0 + T*k) - cos(a0 + k*t4)))/pow(k,2) + (u2*(cos(a0 + k*t1) - cos(a0 + k*t2)))/pow(k,2) + (u1*(cos(a0 + k*t2) - cos(a0 + k*t3)))/pow(k,2) + (u2*(cos(a0 + k*t3) - cos(a0 + k*t4)))/pow(k,2) - (u1*sin(a0 + k*t4)*(T - t4))/k + (u2*sin(a0 + k*t1)*(t1 - t2))/k + (u1*sin(a0 + k*t2)*(t2 - t3))/k + (u2*sin(a0 + k*t3)*(t3 - t4))/k),2);
            double vz_error=lmd4*pow((x4tf - x40 - (u1*(sin(a0 + T*k) - sin(a0 + k*t4)))/k + (u2*(sin(a0 + k*t1) - sin(a0 + k*t2)))/k + (u1*(sin(a0 + k*t2) - sin(a0 + k*t3)))/k + (u2*(sin(a0 + k*t3) - sin(a0 + k*t4)))/k - (u1*(sin(a0 + k*t1) - sin(a0)))/k),2);
            if(y_error+z_error)
            {
                costdata.lmd1=y_error/(y_error+z_error);
                costdata.lmd3=10*z_error/(y_error+z_error);
            }
            if(vy_error+vz_error)
            {
                costdata.lmd2=vy_error/(vy_error+vz_error);
                costdata.lmd4=10*vz_error/(vy_error+vz_error);
            }

//            Eigen::MatrixXf switchMatrix=MatrixXf::Zero(2,array_length+1);// +1 is to avoid the index over in switchMatrix.block(1,array_length,p,q) block operation
//            Eigen::Vector2f inputset(upperthrust,lowthrust);
//            switchMatrix.block(0,0,1,timenodevector(0))=Eigen::MatrixXf::Ones(1,timenodevector(0));
////            if(timenodevector(1)-timenodevector(0))
//            {
//                switchMatrix.block(1,timenodevector(0),1,timenodevector(1)-timenodevector(0))=Eigen::MatrixXf::Ones(1,timenodevector(1)-timenodevector(0));
//            }
//            switchMatrix.block(0,timenodevector(1),1,timenodevector(2)-timenodevector(1))=Eigen::MatrixXf::Ones(1,timenodevector(2)-timenodevector(1));
//            switchMatrix.block(1,timenodevector(2),1,timenodevector(3)-timenodevector(2))=Eigen::MatrixXf::Ones(1,timenodevector(3)-timenodevector(2));
//            switchMatrix.block(0,timenodevector(3),1,array_length-timenodevector(3))=Eigen::MatrixXf::Ones(1,array_length-timenodevector(3));
//            compensatethrust=(inputset * switchMatrix).array().block(0,0,1,array_length);
        }
        ocp_thrust_array=thrust_array+compensatethrust;
        ocp_startindex=0;
//        cout<<"thrust_array:"<<thrust_array<<" compensatethrust: "<<compensatethrust<<" ocp_thrust_array:"<<ocp_thrust_array<<endl;
        return true;
    }
    return false;
}

///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>主 函 数<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
int main(int argc, char **argv)//argc  argument count 传参个数，argument value
{
    ros::init(argc, argv, "position_control");//初始化节点名称
    ros::NodeHandle nh;


    ros::ServiceClient arming_client = nh.serviceClient<mavros_msgs::CommandBool>(
            "mavros/cmd/arming"); //使能解锁飞机  创建client对象并向arming发出请求，服务类型为CommandBool
    ros::ServiceClient setmode_client = nh.serviceClient<mavros_msgs::SetMode>(
            "mavros/set_mode"); //设置为自动控制模式 服务类型为SetMode

    // 【订阅】无人机当前状态/位置/速度信息
    ros::Subscriber state_sub = nh.subscribe<mavros_msgs::State>("mavros/state", 1,
                                                                 state_cb);//订阅器，订阅接收mavros/state话题的mavros_msgs::State类型的消息，有消息到达这个话题时会自动调用state_cb函数
    ros::Subscriber plane_position_pose_sub = nh.subscribe<nav_msgs::Odometry>("mavros/local_position/odom", 1,
                                                                               plane_pos_cb);//pos+twist
    ros::Subscriber plane_poseimu_sub = nh.subscribe<sensor_msgs::Imu>("/mavros/imu/data", 1,
                                                                       plane_imu_cb);//pos+twist
    ros::Subscriber plane_velocity_sub = nh.subscribe<geometry_msgs::TwistStamped>(
            "mavros/local_position/velocity_local", 1, plane_vel_cb); //twist
    ros::Subscriber car_position_sub = nh.subscribe<nav_msgs::Odometry>("odom", 1, car_pos_cb); //车的pos+twist
    ros::Subscriber current_pos_sub=nh.subscribe<nav_msgs::Odometry>("currenttarget_postwist",1,wall_twist_cb);

    ros::Subscriber plane_alt_sub = nh.subscribe<std_msgs::Float64>("mavros/global_position/rel_alt", 1, plane_alt_cb);
    ros::Subscriber controlstate_sub = nh.subscribe<offb_posctl::controlstate>("ocp/control_state", 1, controlstate_cb);
    // 【发布】飞机姿态/拉力信息 坐标系:NED系
    ros::Publisher ocplan_postwist_pub = nh.advertise<nav_msgs::Odometry>("ocplan_positiontwist", 1);//bvp计算的期望位置
    ros::Publisher target_atti_thrust_pub = nh.advertise<mavros_msgs::AttitudeTarget>("mavros/setpoint_raw/attitude",
                                                                                      1);//发布给mavros的控制量(经过换算)
    ros::Publisher plane_rpy_pub = nh.advertise<geometry_msgs::Vector3>("drone/current_rpy", 1);//飞机当前的rpy角
    ros::Publisher current_relativepostwist_pub = nh.advertise<nav_msgs::Odometry>("current_relative_postwist",
                                                                                   1);//当前状态方程中的状态量,即相对量
    ros::Publisher targeterror_pub = nh.advertise<geometry_msgs::Vector3>("targeterror", 1);//目标误差
    ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("atucaltrajectory",1, true);

    int constraintsnumber = 3;
    std::vector<double> lb(timeswitchnumber, 0.0);
    std::vector<double> ub(timeswitchnumber, costdata.T);
    std::vector<double> constraints_tol(constraintsnumber, 1e-6);
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_min_objective(myvfunc, &costdata);
    opt.add_inequality_mconstraint(multiconstraint, NULL, constraints_tol);
    opt.set_xtol_rel(1e-4);
    opt.set_maxtime(0.01);

    // 频率 [30Hz]
    ros::Rate rate(controlfreq);   //50hz的频率发送/接收topic  ros与pixhawk之间,50Hz control frequency

    // log输出文件初始化
    logfile.open("/home/sensenliu/catkin_ws/src/gazebo_ros_learning/offb_posctl/data/savetest.csv", std::ios::out);
    if (!logfile.is_open()) {
        ROS_ERROR("log to file error!");
//        return 0;
    }

    // 读取PID参数
    std::string paraadr("/home/sensenliu/catkin_ws/src/gazebo_ros_learning/offb_posctl/data/param");
    if (param.readParam(paraadr.c_str()) == 0) {
        std::cout << "read config file error!" << std::endl;
//        return 0;
    }


    /// 设置速度环PID参数 比例参数 积分参数 微分参数
    PIDVX.setPID(param.vx_p, param.vx_i, param.vx_d);
    PIDVY.setPID(param.vy_p, param.vy_i, param.vy_d);
    PIDVZ.setPID(param.vz_p, param.vz_i, param.vz_d);
    // 设置速度环积分上限 控制量最大值 误差死区
    PIDVX.set_sat(6, 10, 0);
    PIDVY.set_sat(2, 3, 0);
    PIDVZ.set_sat(2, 5, 0);

    /// 等待和飞控的连接
    while (ros::ok() && current_state.connected == 0) {

        ros::spinOnce(); //调用回调函数
        ros::Duration(1).sleep();
        ROS_INFO("Not Connected");
    }
    ROS_INFO("Connected!!");

    target_atti_thrust_msg.orientation.x = 0;
    target_atti_thrust_msg.orientation.y = 0;
    target_atti_thrust_msg.orientation.z = 0;
    target_atti_thrust_msg.orientation.w = -1;
    target_atti_thrust_msg.thrust = 0.0; //65%的油门 50%与重力平衡即悬停

    /// get car current pose to set plane pose
    float x = pose_car_odom.pose.pose.orientation.x;
    float y = pose_car_odom.pose.pose.orientation.y;
    float z = pose_car_odom.pose.pose.orientation.z;
    float w = pose_car_odom.pose.pose.orientation.w;
    Yaw_Init = quaternion2euler(x, y, z, w).z;

    ///set initial hover position and pose
    target_atti_thrust_msg.orientation.x = x;
    target_atti_thrust_msg.orientation.y = y;
    target_atti_thrust_msg.orientation.z = z;
    target_atti_thrust_msg.orientation.w = w;
    ROS_INFO("got initial point ");

    for (int i = 10; ros::ok() && i >0; --i)// let drone take off slightly at begining, but this step seems useless because the drono has not been armed
    {
        target_atti_thrust_pub.publish(target_atti_thrust_msg);
        ros::spinOnce();//让回调函数有机会被执行
        rate.sleep();
    }
    ROS_INFO("OUT OF LOOP WAIT");


    mavros_msgs::CommandBool arm_cmd; //解锁
    mavros_msgs::SetMode offb_set_mode;
    offb_set_mode.request.custom_mode = "OFFBOARD";
    arm_cmd.request.value = true;

    ros::Time last_request = ros::Time::now();

    ///解锁飞机
    while (ros::ok()) {
        if (current_state.mode != "OFFBOARD" &&
            (ros::Time::now() - last_request > ros::Duration(5.0))) {
            if (setmode_client.call(offb_set_mode) && offb_set_mode.response.mode_sent) {
                ROS_INFO("Offboard enabled");
            }
            last_request = ros::Time::now();
        } else {
            if (!current_state.armed &&
                (ros::Time::now() - last_request > ros::Duration(5.0))) {
                if (arming_client.call(arm_cmd) &&
                    arm_cmd.response.success) {
                    ROS_INFO("Vehicle armed");
                    break;
                }
                last_request = ros::Time::now();
            }
        }
        target_atti_thrust_pub.publish(target_atti_thrust_msg);
        ros::spinOnce();
        rate.sleep();
    }


    /// reach initial hover position and pose by position control
    /*float error_position_sum = 0;
    float error_pose_sum = 5;
    float yaw_current;
    ros::Time begin_time_01 = ros::Time::now();
    int count = 0;
    float thrust_target_sum = 0;
    vector<float> orientation_x;
    vector<float> orientation_y;
    vector<float> orientation_z;
    vector<float> orientation_w;
    got_initial_point = true;
    base_atti_thrust_msg.thrust=0.57;

    while (ros::ok() && count < 500)//没有给小车速度时始终在这个循环里
    {

        ros::spinOnce();
        plane_expected_position.x = 0; //-1 means axis difference
        plane_expected_position.y = 0; //-1 means axis difference
        plane_expected_position.z = 0.5;

        float cur_time_01 = get_ros_time(begin_time_01);  // 相对时间
        pix_controller(cur_time_01);
        target_atti_thrust_msg.header.stamp.sec = pose_car_odom.header.stamp.sec;
        target_atti_thrust_msg.header.stamp.nsec = pose_car_odom.header.stamp.nsec;
        target_atti_thrust_msg.orientation = orientation_target;
        target_atti_thrust_msg.thrust = thrust_target;

        temp_angle = quaternion2euler(pose_drone_odom.pose.pose.orientation.x, pose_drone_odom.pose.pose.orientation.y,
                                      pose_drone_odom.pose.pose.orientation.z, \
                                      pose_drone_odom.pose.pose.orientation.w);//欧拉角
        rpy.y = temp_angle.y;
        rpy.x = temp_angle.x;
        rpy.z = temp_angle.z;
        plane_rpy_pub.publish(rpy);


        target_atti_thrust_pub.publish(target_atti_thrust_msg);
        rate.sleep();//休息

        count += 1;
        if (count >= 500)//count到达500以后,不会再增加
        {
            cout << "You can run Waffle pi!" << endl;
        }

    }

    ROS_INFO("reached initial point and pose ");
    */
    // 记录启控时间
    ros::Time begin_time_02 = ros::Time::now();
    int tempcounter=0;
    float ascentvel=15;
    float tempCurrentPx=0,tempCurrentPy=0,tempCurrentPz=0;
    float tempgoalPx=0,tempgoalPy=0,tempgoalPz=0;
    int lefnodeindex = 0;
    int rightnodeindex = 0;

///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>主  循  环<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    nav_msgs::Path actual_path;
    geometry_msgs::PoseStamped actual_pose_stamped;
    actual_path.header.stamp=ros::Time::now();
    actual_path.header.frame_id="ground_link";
    int pathseq=0;
    double hoverheight=0.6;

    Vector3d dronetowall_pos_inworld(0,0,0);
    Vector3d dronetowall_vel_inworld(0,0,0);
    Vector3d wall_zaxis(0,0,0);
    Vector3d wall_xaxis(1,0,0);
    Vector3d wall_yaxis(0,0,0);
    while (ros::ok()) {
        float cur_time = get_ros_time(begin_time_02);  // 当前时间

        ax_ini= ax_sgFilter.sgfilter(ax_ini);
        ay_ini= ay_sgFilter.sgfilter(ay_ini);
        az_ini= az_sgFilter.sgfilter(az_ini);
        x_acccomp=x_dob.acc_dob_output(ax_ini,ax_des);
        y_acccomp=y_dob.acc_dob_output(ay_ini,ay_des);
        z_acccomp=z_dob.acc_dob_output(az_ini,az_des);
        if (planeupdateflag)   //订阅到飞机位置则flag为true，发布完相对位置的消息后flag置false
        {
//            ay_ini=min(max((double)ay_ini,-2*9.8),2*9.8);
//            az_ini=min(max((double)az_ini,-9.8),9.8);
            if(planstopflag==false && quad_state>0)
//            if(planstopflag==false)
            {
                px_ini = pose_drone_odom.pose.pose.position.x;
                py_ini = pose_drone_odom.pose.pose.position.y;
                pz_ini = pose_drone_odom.pose.pose.position.z;

                vx_ini = vel_drone.twist.linear.x;
                vy_ini = vel_drone.twist.linear.y;
                vz_ini = vel_drone.twist.linear.z;

//                pz_ini=min(max((double)pz_ini, 0.2),2.5);
//                vz_ini=min(max((double)vz_ini, -5.0),5.0);
//                vy_ini=min(max((double)vy_ini, -5.0),5.0);

                current_relativepostwist_msg.pose.pose.position.x = px_ini;
                current_relativepostwist_msg.pose.pose.position.y = py_ini;
                current_relativepostwist_msg.pose.pose.position.z = pz_ini;

                current_relativepostwist_msg.twist.twist.linear.x = vx_ini;
                current_relativepostwist_msg.twist.twist.linear.y = vy_ini;
                current_relativepostwist_msg.twist.twist.linear.z = vz_ini;

                current_relativepostwist_msg.twist.twist.angular.y = ay_ini;
                current_relativepostwist_msg.twist.twist.angular.z = az_ini;
//            current_relativepostwist_msg.pose.pose.orientation.x=phi_ini;
//            current_relativepostwist_msg.pose.pose.orientation.y=omegax_ini;
//            current_relativepostwist_msg.pose.pose.orientation.z=thrust_ini;
//            current_relativepostwist_msg.pose.pose.orientation.w=tau_ini;
                current_relativepostwist_msg.header.stamp = pose_drone_odom.header.stamp;
                current_relativepostwist_pub.publish(current_relativepostwist_msg);
            }

            actual_pose_stamped.pose.position.x=pose_drone_odom.pose.pose.position.x;
            actual_pose_stamped.pose.position.y=pose_drone_odom.pose.pose.position.y;
            actual_pose_stamped.pose.position.z=pose_drone_odom.pose.pose.position.z;
//            actual_pose_stamped.pose.orientation=pose_drone_odom.pose.pose.orientation;
            actual_pose_stamped.header.stamp=ros::Time::now();
            actual_pose_stamped.header.frame_id="ground_link";
//            actual_pose_stamped.header.seq=i;
            actual_path.poses.push_back(actual_pose_stamped);
            path_pub.publish(actual_path);

//            std::cout << "py_ini:  " << py_ini << " pz_ini:  " << pz_ini << " phi_ini:  " << phi_ini << " vy_ini:  "<< vy_ini <<" vz_ini: "<<vz_ini<<std::endl;//输出,放到py文件中求解
            planeupdateflag = false;
        }

        ros::spinOnce();//刷新callback的消息
        ROS_ERROR_STREAM("quad_state"<<quad_state<<endl);
        switch (quad_state) {
            case 0:
                plane_expected_position.z=plane_expected_position.z+ascentvel*cur_time;
                plane_expected_position.z=min(plane_expected_position.z,hoverheight);
                plane_expected_position.x=0;
                plane_expected_position.y=0;
                plane_expected_velocity.x=0;
                plane_expected_velocity.y=0;
                plane_expected_velocity.z=0;
                plane_expected_acceleration.x=0;
                plane_expected_acceleration.y=0;
                plane_expected_acceleration.z=0;

                pix_controller(cur_time);

                if(plane_expected_position.z>=hoverheight)
                {
                    tempcounter++;
                    if(tempcounter>=150)
                    {
                        quad_state=1;
                        tempcounter=0;
                    }
                }
                break;
            case 1:
                cout<<"-------controlcounter:"<<controlcounter<<"  controlstatearray_msg.tfnodenumber: "<<controlstatearray_msg.tfnodenumber<<endl;
                if(contstaterecieveflag)  //订阅到bvp计算的控制量则flag为true,用于起始时刻,还没算出bvp时
                {
                    lefnodeindex = controlcounter;
//                    cout<<"-------lefnodeindex"<<lefnodeindex<<endl;
                    if (lefnodeindex < controlstatearray_msg.arraylength)
                    {
//                        Here,we used lefnodeindex-1, it is because we consider that the pos and vel error is in next state, not in current state
                        plane_expected_position.x=controlstatearray_msg.stateXarray[lefnodeindex];
                        plane_expected_position.y=controlstatearray_msg.stateYarray[lefnodeindex];
                        plane_expected_position.z=controlstatearray_msg.stateZarray[lefnodeindex];


                        plane_expected_velocity.x=controlstatearray_msg.stateVXarray[lefnodeindex];
                        plane_expected_velocity.y=controlstatearray_msg.stateVYarray[lefnodeindex];
                        plane_expected_velocity.z=controlstatearray_msg.stateVZarray[lefnodeindex];

                        plane_expected_acceleration.x=controlstatearray_msg.stateAXarray[lefnodeindex];
                        plane_expected_acceleration.y=controlstatearray_msg.stateAYarray[lefnodeindex];
                        plane_expected_acceleration.z=controlstatearray_msg.stateAZarray[lefnodeindex];
                    }else
                    {
                        plane_expected_position.x=controlstatearray_msg.stateXarray[controlstatearray_msg.stateXarray.size()-1];
                        plane_expected_position.y=controlstatearray_msg.stateYarray[controlstatearray_msg.stateXarray.size()-1];
                        plane_expected_position.z=controlstatearray_msg.stateZarray[controlstatearray_msg.stateXarray.size()-1];


                        plane_expected_velocity.x=controlstatearray_msg.stateVXarray[controlstatearray_msg.stateXarray.size()-1];
                        plane_expected_velocity.y=controlstatearray_msg.stateVYarray[controlstatearray_msg.stateXarray.size()-1];
                        plane_expected_velocity.z=controlstatearray_msg.stateVZarray[controlstatearray_msg.stateXarray.size()-1];

                        plane_expected_acceleration.x=controlstatearray_msg.stateAXarray[controlstatearray_msg.stateXarray.size()-1];
                        plane_expected_acceleration.y=controlstatearray_msg.stateAYarray[controlstatearray_msg.stateXarray.size()-1];
                        plane_expected_acceleration.z=controlstatearray_msg.stateAZarray[controlstatearray_msg.stateXarray.size()-1];
                    }
                    controlcounter++;
                    if((controlstatearray_msg.arraylength-controlcounter)<0.3*controlstatearray_msg.discrepointpersecond && controlstatearray_msg.tfnodenumber==2)
                    {
                        startattitudecotrolflag=true;
                        ROS_ERROR_STREAM("startattitudecotrolflag:"<<startattitudecotrolflag<<" plane_expected_position.y:"<<plane_expected_position.y<<" plane_expected_position.z: "<<plane_expected_position.z<<" current y: "<<pose_drone_odom.pose.pose.position.y<<" current z: "<<pose_drone_odom.pose.pose.position.z<<" planned roll:"<<angle_target.x<<" current roll:"<<angle_receive.x);
                    }
                    if(controlcounter>=controlstatearray_msg.arraylength && controlstatearray_msg.tfnodenumber==2)
                    {
                        ROS_ERROR_STREAM( "end_state2---plane_expected_position.y:"<<plane_expected_position.y<<" plane_expected_position.z:"<<plane_expected_position.z<<" current y: "<<pose_drone_odom.pose.pose.position.y<<" current z: "<<pose_drone_odom.pose.pose.position.z<<" planned roll: "<<-atan(controlstatearray_msg.stateAYarray[lefnodeindex]/(controlstatearray_msg.stateAZarray[lefnodeindex]+9.8))<<" controlstatearray_msg.tfnodenumber: "<<controlstatearray_msg.tfnodenumber
                        <<" controlstatearray_msg.arraylength: "<<controlstatearray_msg.arraylength);
                        quad_state=2;
                        wall_zaxis << 0,controlstatearray_msg.stateAYarray[controlstatearray_msg.arraylength-1]/9.8,(controlstatearray_msg.stateAZarray[controlstatearray_msg.arraylength-1]+9.8)/9.8;
                        wall_yaxis=wall_zaxis.cross(wall_xaxis);

                    }
//                    planned_postwist_msg.pose.pose.orientation.x=-atan(controlstatearray_msg.stateAYarray[lefnodeindex]/(controlstatearray_msg.stateAZarray[lefnodeindex]+9.8));
                }
                pix_controller(cur_time);
                if(startattitudecotrolflag&&targetangleachievedflag== false)
                {
                    terminalreplanning();
                    if(ocp_startindex<ocp_thrust_array.size())
                    {
                        thrust_target=ocp_thrust_array(ocp_startindex)/9.8*(param.THR_HOVER);
                        ocp_startindex++;
                    }
                }
                planned_postwist_msg.pose.pose.position.x=plane_expected_position.x;
                planned_postwist_msg.pose.pose.position.y=plane_expected_position.y;
                planned_postwist_msg.pose.pose.position.z=plane_expected_position.z;
                planned_postwist_msg.twist.twist.linear.y=plane_expected_velocity.y;
                planned_postwist_msg.twist.twist.linear.z=plane_expected_velocity.z;
                planned_postwist_msg.pose.pose.orientation.x=angle_target.x;
                planned_postwist_msg.twist.twist.angular.y=PIDVY.Output;
                planned_postwist_msg.twist.twist.angular.z=PIDVZ.Output;
                ROS_INFO_STREAM("Apporaching_thrust_target: "<< thrust_target<<" roll:"<<angle_target.x<<" pitch:"<<angle_target.y<<" yaw:"<<angle_target.z);
                break;
            case 2:
                tempcounter++;
//                if(tempcounter<=controlstatearray_msg.parabolictime*controlstatearray_msg.discrepointpersecond)
//                {
//                    tempgoalPx=pose_drone_odom.pose.pose.position.x;
//                    tempgoalPy=pose_drone_odom.pose.pose.position.y;
//                    tempgoalPz=pose_drone_odom.pose.pose.position.z;
//                    cout<<"tempcounter: "<<tempcounter<<"  parabolictime*30: "<<controlstatearray_msg.parabolictime*30<<endl;
//                    ROS_ERROR_STREAM("temp goal goal goal Px:"<<tempgoalPx<<" Py: "<<tempgoalPy<<" Pz: "<<tempgoalPz<<" pitch:"<<temp_angle.x);
//                    thrust_target=param.THR_HOVER;
//                    break;
//                }
                tempCurrentPx=pose_drone_odom.pose.pose.position.x;
                tempCurrentPy=pose_drone_odom.pose.pose.position.y;
                tempCurrentPz=pose_drone_odom.pose.pose.position.z;
                if(wallstate_receivedflag)
                {
                    dronetowall_pos_inworld << 0.0, tempCurrentPy-(wallpostwist_odom.pose.pose.position.y),tempCurrentPz-(wallpostwist_odom.pose.pose.position.z);
                    dronetowall_vel_inworld << 0.0, vel_drone.twist.linear.y-wallpostwist_odom.twist.twist.linear.y,vel_drone.twist.linear.z-wallpostwist_odom.twist.twist.linear.z;
                }else{
                    dronetowall_pos_inworld << 0.0, tempCurrentPy-controlstatearray_msg.rendezvouswall_y,tempCurrentPz-controlstatearray_msg.rendezvouswall_z;
                    dronetowall_vel_inworld << 0.0, vel_drone.twist.linear.y-controlstatearray_msg.rendezvouswall_vy,vel_drone.twist.linear.z-controlstatearray_msg.rendezvouswall_vy;

                }
                thrust_target  = param.THR_HOVER*cos(-atan(controlstatearray_msg.stateAYarray[lefnodeindex]/(controlstatearray_msg.stateAZarray[lefnodeindex]+9.8)))*cos(0);   //目标推力值 to alleviate the gravity's component along the drone's z axis
                thrust_target  = max((double)thrust_target,param.THR_HOVER*0.5);   //目标推力值,只是用来保证提供扭矩，the drone is easy to fall freely and crash

                if(targetangleachievedflag== false)
                {
                    terminalreplanning();
                    if(ocp_startindex<ocp_thrust_array.size())
                    {
                        thrust_target=ocp_thrust_array(ocp_startindex)/9.8*(param.THR_HOVER);
                        ocp_startindex++;
                    }
                }
                ROS_ERROR_STREAM(" current y: "<<tempCurrentPy<<" current vy"<<vel_drone.twist.linear.y<<" current z: "<<tempCurrentPz<<" current vz"<<vel_drone.twist.linear.z<<" current roll:"<<angle_receive.x);
                ROS_ERROR_STREAM(" planned rendez y: "<<ytf_ocp<<" ren vy:"<<vytf_ocp<<" ren z:"<<ztf_ocp<<" ren vz:"<<vztf_ocp<<" remaining T:"<<costdata.T);
                ROS_ERROR_STREAM(" current wall y: "<<wallpostwist_odom.pose.pose.position.y<<" current wall vy: "<<wallpostwist_odom.twist.twist.linear.y<<" current wall z:"<<wallpostwist_odom.pose.pose.position.z<<" current wall vz:"<<wallpostwist_odom.twist.twist.linear.z);
                if(dronetowall_pos_inworld.dot(wall_zaxis)<-0.2||dronetowall_vel_inworld.dot(wall_yaxis)<-1.5)
                {
                    quad_state=3;
                    tempcounter=0;
                    ROS_ERROR_STREAM("dronetowall_pos_inworld.dot(wall_zaxis):"<<dronetowall_pos_inworld.dot(wall_zaxis)<<" dronetowall_vel_inworld.dot(wall_yaxis): "<<dronetowall_vel_inworld.dot(wall_yaxis));
                }
                if(tempcounter>=30&&(dronetowall_pos_inworld.dot(wall_zaxis)>=-0.2&&dronetowall_vel_inworld.dot(wall_yaxis)>=-1.5))
//                if(tempcounter>=30&&tempCurrentPz>=controlstatearray_msg.wall_z-0.15)
                {
                    quad_state=4;
                    tempcounter=0;
                    ROS_ERROR_STREAM("Success-tempCurrentPz:"<<tempCurrentPz<<" rendezvouswall_z: "<<controlstatearray_msg.rendezvouswall_z<<" vel_read.z: "<<vel_drone.twist.linear.z);
                }
                orientation_target = euler2quaternion(-atan(controlstatearray_msg.stateAYarray[lefnodeindex]/(controlstatearray_msg.stateAZarray[lefnodeindex]+9.8)), 0, angle_target.z);

                ROS_INFO_STREAM("Duringsuck_thrust_target: "<< thrust_target<<" roll:"<<-atan(controlstatearray_msg.stateAYarray[lefnodeindex]/(controlstatearray_msg.stateAZarray[lefnodeindex]+9.8))<<" pitch:"<<0<<" yaw:"<<angle_target.z);
                break;
            case 3:
                startattitudecotrolflag=false;
                plane_expected_position.x=tempCurrentPx;
                plane_expected_position.y=tempCurrentPy;
                plane_expected_position.z=tempCurrentPz;
                plane_expected_velocity.x=0;
                plane_expected_velocity.y=0;
                plane_expected_velocity.z=0;
                plane_expected_acceleration.x=0;
                plane_expected_acceleration.y=0;
                plane_expected_acceleration.z=0;
                pix_controller(cur_time);
                break;
            case 4:
                orientation_target = euler2quaternion(-atan(controlstatearray_msg.stateAYarray[lefnodeindex]/(controlstatearray_msg.stateAZarray[lefnodeindex]+9.8)), 0, angle_target.z);
                thrust_target  = 0;   //目标推力值
                ROS_INFO_STREAM("Sucksuccess_thrust_target: "<< thrust_target<<" roll:"<<angle_target.x<<" pitch:"<<angle_target.y<<" yaw:"<<angle_target.z);
                break;
            default:
                startattitudecotrolflag=false;
                break;
        }

        cout<<"refpos x: "<<plane_expected_position.x<<"   y:"<<plane_expected_position.y<<"  z:"<<plane_expected_position.z<<endl;

            ///publish plane current rpy
            temp_angle = quaternion2euler(pose_drone_odom.pose.pose.orientation.x,
                                          pose_drone_odom.pose.pose.orientation.y,
                                          pose_drone_odom.pose.pose.orientation.z,
                                          pose_drone_odom.pose.pose.orientation.w);//欧拉角
            plane_rpy_pub.publish(temp_angle);
//            cout<<"drone_roll--"<<temp_angle.x<<endl;

            ///publish thrust & orientation
//            std::cout << "thrust_target: " << thrust_target << std::endl;
            target_atti_thrust_msg.header.stamp = ros::Time::now();
            target_atti_thrust_msg.orientation = orientation_target;
            target_atti_thrust_msg.thrust = thrust_target;
            target_atti_thrust_pub.publish(target_atti_thrust_msg);

            ///publish planned pos&twist in x&z
            planned_postwist_msg.header.stamp = ros::Time::now();
            ocplan_postwist_pub.publish(planned_postwist_msg);



            ///publish targeterror_msg
            targeterror_msg.x = px_ini + 3;
            targeterror_msg.z = pz_ini + px_ini * rpy.y;
            targeterror_msg.y = py_ini;//pub time consumption
            targeterror_pub.publish(targeterror_msg);
            if(temp_angle.x-1.22>=-0.2)
            {
                anglereachedflag= true;
            }

            if(param.Enable_log_to_file)
            {
                data_log(logfile, cur_time);                     //log输出
            }
            rate.sleep();

        }
        logfile.close();
        return 0;
    }

/**
 * 获取当前时间 单位：秒
 */
float get_ros_time(ros::Time time_begin)
{
    ros::Time time_now = ros::Time::now();
    float currTimeSec = time_now.sec-time_begin.sec;
    float currTimenSec = time_now.nsec / 1e9 - time_begin.nsec / 1e9;
    return (currTimeSec + currTimenSec);
}

///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>函 数 定 义<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
float pitch_compensation(float theta_dc, float theta_dn, float yk[],float uk[])
{
    float y_next,pitch_revise,current_error;

    y_next = 1.461 * yk[0] - 0.5345 * yk[1] + 0.01827 * uk[3] +0.03653 * uk[4] + 0.01827 * uk[5]; // the predicted pitch in the next moment
    current_error = theta_dn - y_next;//current error
    comp_integrate += current_error;//error sum

    pitch_revise = theta_dc + comp_kp*current_error + comp_ki*comp_integrate + comp_kd*(current_error - comp_last_error)*controlfreq;
    comp_last_error = current_error;

    return (pitch_revise);
}

void vector3dLimit(Vector3d &v, double limit)  ///limit should be positive
{
    if(limit > 0){
        for(int i=0; i<3; i++){
            v(i) = fabs(v(i)) > limit ? (v(i) > 0 ? limit : -limit) : v(i);
        }
    }
}

Vector3d vectorElementMultiply(Vector3d v1, Vector3d v2)
{
    Vector3d result;
    result << v1(0)*v2(0), v1(1)*v2(1), v1(2)*v2(2);
    return result;
}

int pix_controller(float cur_time)
{
//位 置 环
    //计算误差
    float error_x = plane_expected_position.x - pose_drone_odom.pose.pose.position.x;
    float error_y = plane_expected_position.y - pose_drone_odom.pose.pose.position.y;
    float error_z = plane_expected_position.z - plane_real_alt.data;
//    std::cout << "error: x：" << error_x << "\ty：" << error_y << "\tz：" << error_z << std::endl;
    float acc_xd = param.x_p * error_x;
    float acc_yd = param.y_p * error_y;
    float acc_zd = param.z_p * error_z;


    //计算误差
    float error_vx = plane_expected_velocity.x - vel_drone.twist.linear.x;
    float error_vy = plane_expected_velocity.y  - vel_drone.twist.linear.y;
    float error_vz = plane_expected_velocity.z  - vel_drone.twist.linear.z;
    float acc_vxd = param.vx_p * error_vx;
    float acc_vyd = param.vy_p * error_vy;
    float acc_vzd = param.vz_p * error_vz;

    feedforwardcoefx=1/(1+fabs(acc_xd+acc_vxd)/(fabs(plane_expected_acceleration.x)+0.5));//0.5 is to avoid the zero of denominator
    feedforwardcoefy=1/(1+fabs(acc_yd+acc_vyd)/(fabs(plane_expected_acceleration.y)+0.5));
    feedforwardcoefz=1/(1+fabs(acc_zd+acc_vzd)/(fabs(plane_expected_acceleration.z)+0.5));
//    cout<<"feedforwardcoefx:"<<feedforwardcoefx<<" y:"<<feedforwardcoefy<<" z:"<<feedforwardcoefz<<endl;
    //计算输出
    if(startattitudecotrolflag==false)
    {
        PIDVX.Output=acc_xd+acc_vxd+feedforwardcoefx*plane_expected_acceleration.x+x_acccomp;
        PIDVY.Output=acc_yd+acc_vyd+feedforwardcoefy*plane_expected_acceleration.y+y_acccomp;
        PIDVZ.Output=acc_zd+acc_vzd+feedforwardcoefz*plane_expected_acceleration.z+z_acccomp;
    }else{
        PIDVX.Output=plane_expected_acceleration.x;
        PIDVY.Output=plane_expected_acceleration.y;
        PIDVZ.Output=plane_expected_acceleration.z;
    }
//    cout<<"acc_p_error  planz_a:"<<plane_expected_acceleration.z<<"  vz_a:"<<acc_vzd<<"  z:"<<acc_zd<<" roll: "<<angle_receive.x<<" PIDVZ.Output: "<<PIDVZ.Output<<endl;

    angle_target.x = asin(-PIDVY.Output/sqrt(pow(PIDVX.Output,2)+pow(PIDVY.Output,2)+pow(PIDVZ.Output+9.8,2)));
    angle_target.y = atan(PIDVX.Output/(PIDVZ.Output+9.8));
    angle_target.z = Yaw_Init;

    orientation_target = euler2quaternion(angle_target.x, angle_target.y, angle_target.z);
//    thrust_target = (float)(0.05 * (9.8 + PIDVZ.Output));   //目标推力值
    thrust_target  = (float)sqrt(pow(PIDVX.Output,2)+pow(PIDVY.Output,2)+pow(PIDVZ.Output+9.8,2))/9.8*(param.THR_HOVER);   //目标推力值

//    std::cout << "PIDVZ.OUTPUT:  " << PIDVZ.Output << std::endl;
//    std::cout << "thrust_target:  " << thrust_target <<"param.THR_HOVER:"<<param.THR_HOVER<<std::endl;

    return 0;
}
/**
 * 将欧拉角转化为四元数
 * @param roll
 * @param pitch
 * @param yaw
 * @return 返回四元数
 */
geometry_msgs::Quaternion euler2quaternion(float roll, float pitch, float yaw)
{
    geometry_msgs::Quaternion temp;
    temp.w = cos(roll/2)*cos(pitch/2)*cos(yaw/2) + sin(roll/2)*sin(pitch/2)*sin(yaw/2);
    temp.x = sin(roll/2)*cos(pitch/2)*cos(yaw/2) - cos(roll/2)*sin(pitch/2)*sin(yaw/2);
    temp.y = cos(roll/2)*sin(pitch/2)*cos(yaw/2) + sin(roll/2)*cos(pitch/2)*sin(yaw/2);
    temp.z = cos(roll/2)*cos(pitch/2)*sin(yaw/2) - sin(roll/2)*sin(pitch/2)*cos(yaw/2);
    return temp;
}

/**
 * 将四元数转化为欧拉角形式
 * @param x
 * @param y
 * @param z
 * @param w
 * @return 返回Vector3的欧拉角
 */
geometry_msgs::Vector3 quaternion2euler(float x, float y, float z, float w)
{
    geometry_msgs::Vector3 temp;
    temp.x = atan2(2.0 * (w * x + y * z), 1.0 - 2.0 * (x * x + y * y));
    temp.y = asin(2.0 * (w * y - z * x));
    temp.z = atan2(2.0 * (w * z + x * y), 1.0 - 2.0 * (y * y + z * z));
    return temp;
}

/**
 * 将进入offboard后的位置&速度&姿态信息记录进文件
 * @param cur_time
 */

void data_log(std::ofstream &logfile, float cur_time)
{
    logfile << cur_time << ","
            <<plane_expected_position.x <<","<<plane_expected_position.y <<","<<plane_expected_position.z <<","    //set_pos
            <<pos_drone.pose.position.x <<","<<pos_drone.pose.position.y <<","<<pos_drone.pose.position.z <<","    //uav_pos
            <<plane_expected_velocity.x <<","<<plane_expected_velocity.y <<","<<plane_expected_velocity.z <<","                                           //set_vel
            <<vel_read.x <<","<<vel_read.y <<","<<vel_read.z <<","       //uav_vel
            <<PIDVX.Output <<","<<PIDVY.Output<<","<<PIDVZ.Output <<","       //uav_vel
            <<0 <<","<<ay_ini<<","<<az_ini <<","       //uav_vel
            <<angle_target.x  <<","<<angle_target.y  <<","<<angle_target.z  <<","                                  //set_att
            <<angle_receive.x <<","<<angle_receive.y <<","<<angle_receive.z <<","
            <<pose_drone_Imu.angular_velocity.x <<","<<pose_drone_Imu.angular_velocity.y <<","<<pose_drone_Imu.angular_velocity.z <<","
            <<pos_wall.pose.position.y<<","<<pos_wall.pose.position.z<<","<<vel_wall.twist.linear.y<<","<<vel_wall.twist.linear.z<<","
            <<controlstatearray_msg.predictstartwall_y<<","<<controlstatearray_msg.predictstartwall_z<<","<<controlstatearray_msg.predictstartwall_vy<<","<<controlstatearray_msg.predictstartwall_vz<<","//wall prediction start point
            <<controlstatearray_msg.rendezvouswall_y<<","<<controlstatearray_msg.rendezvouswall_z<<","<<controlstatearray_msg.rendezvouswall_vy<<","<<controlstatearray_msg.rendezvouswall_vz<<","//wall prediction rendezvous point
            <<controlstatearray_msg.timecost<<","<<controlstatearray_msg.arraylength<<","<<controlstatearray_msg.increaseratio<<","//calcualte time, path time
            <<startattitudecotrolflag<<","<<anglereachedflag<<","<<feedforwardcoefx<<","<<feedforwardcoefy<<","<<feedforwardcoefz<<","//calcualte time, path time
            <<thrust_target<<std::endl;
//                cout<<"pos_wall.pose.position.y: "<<pos_wall.pose.position.y<<" pos_wall.header.stamp:"<< fixed << setprecision(5)<<pos_wall.header.stamp.toSec()<<" controlstatearray_msg.header.stamp:"<< fixed << setprecision(5)<< controlstatearray_msg.header.stamp.toSec()<<endl;

}