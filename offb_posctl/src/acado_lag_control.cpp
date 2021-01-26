//
// Created by lynn on 2020/7/19.
//
//#include <acado/acado_optimal_control.hpp>
#include <acado/acado_toolkit.hpp>
#include "ros/ros.h"
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
int pointnumber=20;// the number is almost always 20. It less, the accuracy won't be enough, if more, the time consumpiton will be too large.
int controlfreq=50;
int discretizedpointpersecond=(int)pointnumber/t_end;
float px_ini = -3.0;
float pz_ini = 0.0;
float vx_ini = -0.1;
float vz_ini = 0.0;
float theta_ini = 0.0;
float rate_ini = 0.0;
float vxPlane_ini=0.0,vyPlane_ini=0.0,vzPlane_ini=0.0;
float aircoeffx=1.5,aircoeffz=0.35;
float carPlaneDistance=3.0;
offb_posctl::controlstate controlstate_msg;
bool currentupdateflag= false;
FILTER filterDroneVelx(150);
sensor_msgs::Imu drone_imu;
void do_process()
{
    chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now();
    DifferentialState        px,pz,vx,vz,vxPlane,theta,rate;     // px=xdrone-xcar, pz=zdrone-zcar vx=dot(px),vz=dot(pz)
    Control                  u1,u2          ;     // u1 is thrust, u2 is theta
    DifferentialEquation     f( 0.0, t_end );     // the differential equation
    f << dot(px) == vx;                         // an implementation
    f << dot(pz) == vz;             // of the model equations
    f << dot(vx) == u1*sin(theta)-aircoeffx*vxPlane;                 // for the drone.
    f << dot(vz) == u1*cos(theta)-9.8-aircoeffz*vxPlane;
    f << dot(vxPlane) == u1*sin(theta)-aircoeffx*vxPlane;
    f << dot(theta)==rate;
    f << dot(rate)==u2;
    OCP ocp_(0.0,t_end,pointnumber);
    ocp_.minimizeLagrangeTerm(1*(pz+px*theta)*(pz+px*theta)+3*(px+carPlaneDistance)*(px+carPlaneDistance)+0.5*(u1-9.8)*(u1-9.8)+0.5*u2*u2); // the 1.5 is for the target high(in this gazebo the height of car is 0)
//        ocp_.minimizeLagrangeTerm(1*(pz+px*u2)*(pz+px*u2)+3*(px+carPlaneDistance)*(px+carPlaneDistance)+0.5*(u1-9.8)*(u1-9.8)+0.5*u2*u2); // the 1.5 is for the target high(in this gazebo the height of car is 0)
//    ocp_.minimizeLagrangeTerm((pz-0.61-(0.1*t-px)*u2)*(pz-0.61-(0.1*t-px)*u2)+(0.1*t-px-3)*(0.1*t-px-3)+0.5*u1*u1+0.5*u2*u2);

    ocp_.subjectTo( f                   );     // minimize T s.t. the model,
    ocp_.subjectTo( AT_START, px == px_ini );     // the initial values for s,
    ocp_.subjectTo( AT_START, pz == pz_ini);     // v,
    ocp_.subjectTo( AT_START, vx == vx_ini);     // and m,
    ocp_.subjectTo( AT_START, vz == vz_ini);     // and m,
    ocp_.subjectTo( AT_START, vxPlane == vxPlane_ini);     // and m,
    ocp_.subjectTo( AT_START, theta == theta_ini);     // and m,
    ocp_.subjectTo( AT_START, rate == rate_ini);     // and m,

    ocp_.subjectTo( 0 <= u1 <=  9.8*2   );     // the crol input u,
    ocp_.subjectTo(  -0.628<= theta <= 0.628  );     // and the time horizon T.
    ocp_.subjectTo(  -1<= vxPlane <= 10  );     // and the time horizon T.
    ocp_.subjectTo(  -3<= rate <= 3  );
    ocp_.subjectTo(  -20<= u2 <= 20  );

    OptimizationAlgorithm algorithm(ocp_);     // the optimization algorithm
    algorithm.set(MAX_NUM_ITERATIONS,50);// default value is 1000
    algorithm.set(MAX_NUM_QP_ITERATIONS,100);//default value is 10000
    algorithm.set(KKT_TOLERANCE,1e-3);// default value is 1e-6
    algorithm.set(PRINTLEVEL,NONE);// do not print solution process.
    try {
        resultofsolve=algorithm.solve();                        // solves the problem.
        algorithm.getControls(control);
        algorithm.getDifferentialStates(state);
//        cout << "controls:  time   |  controls" <<endl;
//        control.print();
    }catch (exception e)
    {
        resultofsolve=-1;
        std::cout<<"ocp throw exception"<<std::endl;
    }
    px.clearStaticCounters();
    pz.clearStaticCounters();
    vx.clearStaticCounters();
    vz.clearStaticCounters();
    vxPlane.clearStaticCounters();
    theta.clearStaticCounters();
    rate.clearStaticCounters();
    u1.clearStaticCounters();
    u2.clearStaticCounters();
    chrono::time_point<chrono::steady_clock> end_time = chrono::steady_clock::now();
    timeconsumption=chrono::duration_cast<chrono::milliseconds>(end_time - begin_time).count();
    cout <<"time acado consume~~...........~~~~~~~~mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm: "<<timeconsumption<< endl;
//    controlcounter=floor(timeconsumption/1000.0*controlfreq)
}

/**match the current state with the newest states solved by OCP, and the time of the nearest state will be the first used point
// in the newest solved input. the return value is the matched index
*/
double tempStateX=0.0,tempStateZ=0.0,tempStateVX=0.0,tempStateVZ=0.0;
int StateMatch()
{
    int k1=1,k2=1,k3=1,k4=1,lefnodeindex=0,rightnodeindex=0,controlcounter=0,mincounter=0;
    float minStateDistance=0.0,currentStateDistance=0.0,stateX=0.0,stateZ=0.0,stateVX=0.0,stateVZ=0.0;
    for(;rightnodeindex<controlstate_msg.stateXarray.size();)
    {
        lefnodeindex=floor(controlcounter*discretizedpointpersecond/controlfreq);
        rightnodeindex=lefnodeindex+1;
        if(rightnodeindex<controlstate_msg.stateXarray.size())
        {
            //interpolition
            stateX= (rightnodeindex-(double)controlcounter*discretizedpointpersecond/controlfreq)*controlstate_msg.stateXarray[lefnodeindex]
                                +((double)controlcounter*discretizedpointpersecond/controlfreq-lefnodeindex)*controlstate_msg.stateXarray[rightnodeindex];

            stateZ =(rightnodeindex-(double)controlcounter*discretizedpointpersecond/controlfreq)*controlstate_msg.stateZarray[lefnodeindex]
                               +((double)controlcounter*discretizedpointpersecond/controlfreq-lefnodeindex)*controlstate_msg.stateZarray[rightnodeindex];

            stateVX= (rightnodeindex-(double)controlcounter*discretizedpointpersecond/controlfreq)*controlstate_msg.stateVXarray[lefnodeindex]
                    +((double)controlcounter*discretizedpointpersecond/controlfreq-lefnodeindex)*controlstate_msg.stateVXarray[rightnodeindex];

            stateVZ =(rightnodeindex-(double)controlcounter*discretizedpointpersecond/controlfreq)*controlstate_msg.stateVZarray[lefnodeindex]
                    +((double)controlcounter*discretizedpointpersecond/controlfreq-lefnodeindex)*controlstate_msg.stateVZarray[rightnodeindex];
            currentStateDistance=k1*pow(stateX-px_ini,2)+k2*pow(stateZ-pz_ini,2)+k3*pow(stateVX-vx_ini,2)+k4*pow(stateVZ-vx_ini,2);
            if(minStateDistance==0)
            {
                minStateDistance=currentStateDistance;
                controlcounter=0;
            } else{
                if(minStateDistance>currentStateDistance)
                {
                    minStateDistance=currentStateDistance;
                    mincounter=controlcounter;
                    tempStateX=stateX;
                    tempStateZ=stateZ;
                    tempStateVX=stateVX;
                    tempStateVZ=stateVZ;
                }
            }
        }
        controlcounter++;
    }
//    std::cout<<"matchedcoutner:---"<<mincounter<<" X:"<<tempStateX<<" Z:"<<tempStateZ<<" VX:"<<tempStateVX<<" VZ:"<<tempStateVZ<<endl;
    return mincounter;
}


void state_cb(const mavros_msgs::State::ConstPtr &msg){
    current_state = *msg;
}

bool hasGotImu = false;
void droneImu_cb(const sensor_msgs::Imu::ConstPtr &msg)
{
    drone_imu=*msg;
    rate_ini=drone_imu.angular_velocity.y;
}

void dronerpy_cb(const geometry_msgs::Vector3::ConstPtr &msg)
{
    currentRPYangle =*msg;
    theta_ini=currentRPYangle.y;
//    std::cout <<"currentRPYangle-----------theta_ini:"<<theta_ini<<endl;
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

void plane_vel_cb(const geometry_msgs::TwistStamped::ConstPtr &msg){
    vel_drone = *msg;
//    vxPlane_ini=filterDroneVelx.filter(vel_drone.twist.linear.x);
    vxPlane_ini=(vel_drone.twist.linear.x);
    filteredPlaneVelmsg.x=vxPlane_ini;
    vyPlane_ini=vel_drone.twist.linear.y;
    vzPlane_ini=vel_drone.twist.linear.z;
}

void relative_postwist_cb(const nav_msgs::Odometry::ConstPtr &msg)
{
    px_ini=(*msg).pose.pose.position.x;
    pz_ini=(*msg).pose.pose.position.z;
    vx_ini=(*msg).twist.twist.linear.x;
    vz_ini=(*msg).twist.twist.linear.z;
    currentupdateflag=true;
//    std::cout <<"ocpsolveriniconditons-----------px_ini:  " << px_ini <<"   vx_ini:  " << vx_ini<<"  pz_ini:  " << pz_ini <<"   vz_ini:  " << vz_ini<<
//    "   theta_ini:  " << theta_ini<< "   rate_ini:  " << rate_ini<<std::endl;

}
int main( int argc, char ** argv)
{
    ros::init(argc, argv, "acado_lag_control");
    ros::NodeHandle nh;
    ros::Rate rate(100);// it is the rate to check and process callback functions. the high rate is to avoid the latency to accept the pos and vel
//    //Gazebo 仿真数据
    ros::Subscriber drone_imu_sub=nh.subscribe<sensor_msgs::Imu>("/mavros/imu/data",10,droneImu_cb);
    ros::Subscriber plane_rpy_sub = nh.subscribe<geometry_msgs::Vector3>("drone/current_rpy",10,dronerpy_cb);
    ros::Subscriber car_position_sub = nh.subscribe<nav_msgs::Odometry>("current_relative_postwist",10,relative_postwist_cb); //车的pos+twist
    ros::Subscriber plane_velocity_sub = nh.subscribe<geometry_msgs::TwistStamped>("mavros/local_position/velocity_local", 10, plane_vel_cb); //twist

    // 【发布】飞机姿态/拉力信息 坐标系:NED系
    ros::Publisher controlstate_pub=nh.advertise<offb_posctl::controlstate>("ocp/control_state",10);
    ros::Publisher planeVel_pub=nh.advertise<geometry_msgs::Vector3>("filteredPlaneVel",10);
//  -------------------------------------
    while (ros::ok())
    {
        ros::spinOnce();// to examine the queues of the callback functions once
        if(currentupdateflag)
        {
            do_process();
            if(resultofsolve==SUCCESSFUL_RETURN)
            {
                controlstate_msg.discrepointpersecond=discretizedpointpersecond;
                controlstate_msg.inicounter=0;
                controlstate_msg.arraylength=control.getNumPoints();
                controlstate_msg.thrustarray.clear();
                controlstate_msg.thetaarray.clear();
                controlstate_msg.stateXarray.clear();
                controlstate_msg.stateZarray.clear();
                controlstate_msg.stateVXarray.clear();
                controlstate_msg.stateVZarray.clear();
                for(int i=0;i<control.getNumPoints();i++)
                {
//                    cout<<"control.getNumPoints()-----!!!!:"<<control.getNumPoints()<<"  control.getMatrix(i)(0,0):"<<control.getMatrix(i)(0,0)<<endl;
                    controlstate_msg.thrustarray.push_back((float)control.getMatrix(i)(0,0));
//                    controlstate_msg.thetaarray.push_back((float)control.getMatrix(i)(1,0));
                    controlstate_msg.stateXarray.push_back((float)state.getMatrix(i)(0,0));
                    controlstate_msg.stateZarray.push_back((float)state.getMatrix(i)(1,0));
                    controlstate_msg.stateVXarray.push_back((float)state.getMatrix(i)(2,0));
                    controlstate_msg.stateVZarray.push_back((float)state.getMatrix(i)(3,0));
                    controlstate_msg.thetaarray.push_back((float)state.getMatrix(i)(5,0));

                }
//                ros::spinOnce();// update to newest state to match the state.
//                controlstate_msg.inicounter=StateMatch();//match the state
                controlstate_pub.publish(controlstate_msg);
//                cout<<"controlstate_msg.stateXarray[0]-----------fffffff:"<<controlstate_msg.stateXarray[0]<<"  state.getMatrix(i)(0,0):"<<state.getMatrix(0)(0,0)<<endl;
            }
            currentupdateflag= false;
            planeVel_pub.publish(filteredPlaneVelmsg);
        }
        rate.sleep();
    }

    return 0;
}
