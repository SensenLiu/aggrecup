//
// Created by lynn on 2020/7/19.
//
//#include <acado/acado_optimal_control.hpp>
#include <acado/acado_toolkit.hpp>
#include "ros/ros.h"
#include "std_msgs/Float32.h"
#include <chrono>
#include "iomanip"
#include <acado/acado_gnuplot.hpp>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/Vector3.h>
#include <geometry_msgs/Quaternion.h>
#include <geometry_msgs/Pose.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <nav_msgs/Odometry.h>
#include <mavros_msgs/State.h>
#include <mavros_msgs/AttitudeTarget.h>
#include <sensor_msgs/Imu.h>
#include <std_msgs/Bool.h>
#include <std_msgs/Float32.h>
#include "offb_posctl/controlstate.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include "FILTER.h"
//#include <Eigen/Eigen>// This step will generate error about eigen ,the reason is not clear by far
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <thread>
#include <math.h>
#include <stdlib.h>
#include <FILTER.h>

using namespace Eigen;
using namespace std;

geometry_msgs::PoseStamped pos_ref;         //无人机参考位置

mavros_msgs::State current_state;           //无人机当前状态(mode arm)
sensor_msgs::Imu   imu_drone;               //读入的无人机的IMU信息 包括姿态角和线加速度

geometry_msgs::PoseStamped  pos_drone;      //读入的无人机当前位置
geometry_msgs::PoseStamped  pos_drone_last; //读入的无人机上一次位置
geometry_msgs::PoseStamped  fused_drone; //融合了vicon数据的无人机位姿，注意是NED坐标下的，需要翻转

geometry_msgs::TwistStamped vel_drone;      //读入的无人机当前速度

geometry_msgs::Vector3 acc_receive;         //读入的无人机线加速度
geometry_msgs::Vector3 angle_receive;       //读入的无人机姿态（欧拉角）
geometry_msgs::Vector3 angle_fromvicon;       //读入的vicon计算的sgement姿态（欧拉角）
geometry_msgs::Vector3 angle_fromvicon_qua;       //读入的vicon计算的sgement姿态（欧拉角）
geometry_msgs::Vector3 angle_fromviconinit;  //读入的vicon计算的sgement initial 姿态（欧拉角）
geometry_msgs::Vector3 fusedangle_receive;

geometry_msgs::Quaternion orientation_target;   //发给无人机的姿态指令

geometry_msgs::Vector3 angle_des;            //线性模型输出的理想值
//geometry_msgs::Vector3 angle_dis;            //DOB控制器估计的扰动值
geometry_msgs::Vector3 angle_target;            //经DOB控制器作用后的实际系统输入值
geometry_msgs::Vector3 vel_target, vel_read, vel_read2;
geometry_msgs::Vector3 currentRPYangle;   //欧拉角
geometry_msgs::Vector3 filteredPlaneVelmsg;   //欧拉角
// debug data
geometry_msgs::Vector3 vel_vicon;
const float MAX_POSITION_MEASURE_ERROR = 1;
Eigen::Quaterniond current_quaternion(0,0,0,0);

bool pose_initialized = false;
bool vel_initialized = false;
int posvel_updatesynflag= false;// the variable to make sure the pos and vel used in ocp is synchronized, In reality, we can put pos and vel in one topic to avoid this problem
USING_NAMESPACE_ACADO
VariablesGrid state, parameter, control;// used for output result
returnValue resultofsolve;
int timeconsumption=0;
float t_end=2;
int pointnumber=60;// the number is almost always 20. It less, the accuracy won't be enough, if more, the time consumpiton will be too large.
int controlfreq=30;
int discretizedpointpersecond=(int)pointnumber/t_end;
offb_posctl::controlstate controlstate_msg;
bool currentupdateflag= false;
FILTER filterDroneVelx(150);
sensor_msgs::Imu drone_imu;


void state_cb(const mavros_msgs::State::ConstPtr &msg){
    current_state = *msg;
}


void pos_cb(const geometry_msgs::PoseStamped::ConstPtr &msg){ //由ConstPtr可以看到msg是指针的引用，因此*msg相当于对被引用的指针对象取值
    if (pose_initialized== false)
    {
        pos_drone_last=*msg;
        pose_initialized=true;
        return;
    }
    if(fabs((*msg).pose.position.x - pos_drone_last.pose.position.x) < MAX_POSITION_MEASURE_ERROR &&
       fabs((*msg).pose.position.y - pos_drone_last.pose.position.y) < MAX_POSITION_MEASURE_ERROR &&
       fabs((*msg).pose.position.z - pos_drone_last.pose.position.z) < MAX_POSITION_MEASURE_ERROR)
    {
        pos_drone_last = pos_drone;
        pos_drone = *msg;
        current_quaternion=Eigen::Quaterniond(pos_drone.pose.orientation.w,pos_drone.pose.orientation.x,pos_drone.pose.orientation.y,pos_drone.pose.orientation.z);
       if(pose_initialized&&vel_initialized)
       {
           if (fabs(pos_drone.header.stamp.nsec-vel_drone.header.stamp.nsec)<30e6)
           {
               posvel_updatesynflag=true;
           }
       }
    } else
    {
        pos_drone_last=*msg;
    }
}

int main( int argc, char ** argv) {
    ros::init(argc, argv, "acado_lag_control");
    ros::NodeHandle nh;
    ros::Rate rate(
            100);// it is the rate to check and process callback functions. the high rate is to avoid the latency to accept the pos and vel

    ros::Publisher controlstate_pub = nh.advertise<offb_posctl::controlstate>("ocp/control_state", 10);

/*******fscanf method can not work very well
    char start, end;
    float pz, py, vz, vy, phi, phirate, thrust, torque;
    FILE *fin_state,*fin_u;// the instance can be used to read the file
    fin_state=fopen("/home/sensenliu/catkin_ws/src/gazebo_ros_learning/offb_posctl/data/state.csv","r");// the directory of th data.csv
    while (!feof(fin_state)) {// juddge whether the fin reach the end line of the csv file
        fscanf(fin_state,"%c,%f,%f,%f,%f,%f,%f,%f,%f,%c",&start,&pz,&py,&vz,&vy,&phi, &phirate, &thrust, &torque, &end);
//      fin_state >> line;// read one line data to "line"
        cout<<"every row:"<<pz<<endl;
        counter++; // counter plus plus
    }
    fclose(fin_state);//close！！
    ***********************/

    controlstate_msg.discrepointpersecond=discretizedpointpersecond;
    controlstate_msg.inicounter=0;
    controlstate_msg.arraylength=control.getNumPoints();
    controlstate_msg.thrustarray.clear();
    controlstate_msg.thetaarray.clear();
    controlstate_msg.stateXarray.clear();
    controlstate_msg.stateZarray.clear();
    controlstate_msg.stateVXarray.clear();
    controlstate_msg.stateVZarray.clear();
    float pz, py, vz, vy, phi, phirate, thrust, torque;
    ifstream fin_state, fin_u;// the instance can be used to read the file
    string line,number; //used for store the one line data from csv temporarily
    fin_state.open("/home/sensenliu/catkin_ws/src/gazebo_ros_learning/offb_posctl/data/state.csv");// the directory of th data.csv

    while (getline(fin_state, line)) {
        istringstream readstr(line);
//        cout << "every row:" << line << endl;

        getline(readstr,number,'\t');//here we use '\t' as the seperator, because we find the '\t' in line during debug.
        getline(readstr,number,'\t');//the second col is tiemstamp,

        getline(readstr,number,'\t');
        pz=atof(number.c_str());

        getline(readstr,number,'\t');
        py=atof(number.c_str());

        getline(readstr,number,'\t');
        vz=atof(number.c_str());

        getline(readstr,number,'\t');
        vy=atof(number.c_str());

        getline(readstr,number,'\t');
        phi=atof(number.c_str());


        getline(readstr,number,'\t');// phirate is useless

        getline(readstr,number,'\t');
        thrust=atof(number.c_str());

        getline(readstr,number,'\t');
        torque=atof(number.c_str());

        controlstate_msg.thrustarray.push_back(thrust);
        controlstate_msg.stateXarray.push_back(py);
        controlstate_msg.stateZarray.push_back(pz);
        controlstate_msg.stateVXarray.push_back(vy);
        controlstate_msg.stateVZarray.push_back(vz);
        controlstate_msg.thetaarray.push_back(phi);
    }


//                ros::spinOnce();// update to newest state to match the state.
//                controlstate_msg.inicounter=StateMatch();//match the state
    controlstate_pub.publish(controlstate_msg);
//                cout<<"controlstate_msg.stateXarray[0]-----------fffffff:"<<controlstate_msg.stateXarray[0]<<"  state.getMatrix(i)(0,0):"<<state.getMatrix(0)(0,0)<<endl;

    return 0;
}
