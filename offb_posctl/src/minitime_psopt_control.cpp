//
// Created by lss on 2021/1/23.
//
#include "psopt.h"
#include "ros/ros.h"
#include <chrono>
#include "iomanip"
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

int timeconsumption=0;
float t_end=2;
int pointnumber=10;// the number is almost always 20. It less, the accuracy won't be enough, if more, the time consumpiton will be too large.
int controlfreq=50;

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

Alg  algorithm;
Sol  solution;
Prob problem;
double py0=0,pz0=0.5,phi0=0,vy0=0,vz0=0,omegax0=0,thrust0=9.8,tau0=0;
double pyf=2,pzf=2,phif=1.57,vyf=0.2*sin(phif),vzf=-0.2*cos(phif),omegaxf=0.0,thrustf=9.8*cos(phif),tauf=0;
////x1=y,x2=z,x3=phi,x4=vy,x5=vz,x6=omega,x7=thrust,x8=tau//

MatrixXd x  ;
MatrixXd u  ;
MatrixXd t  ;
MatrixXd x_refined  = MatrixXd::Zero(8,500);

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return tf;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    adouble x7 = states[ 6 ];
    adouble x8 = states[ 7 ];
//    adouble udf = controls[ 0 ];
//    adouble udtau = controls[ 1 ];
    return  0.5*x7*x7+0.5*x8*x8;
//    return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    adouble x1dot, x2dot, x3dot,x4dot, x5dot, x6dot,x7dot, x8dot;

    adouble x1 = states[ 0 ];
    adouble x2 = states[ 1 ];
    adouble x3 = states[ 2 ];
    adouble x4 = states[ 3 ];
    adouble x5 = states[ 4 ];
    adouble x6 = states[ 5 ];
    adouble x7 = states[ 6 ];
    adouble x8 = states[ 7 ];


    adouble udf = controls[ 0 ];
    adouble udtau = controls[ 1 ];

    x1dot = x4;
    x2dot = x5;
    x3dot = x6;
    x4dot = -x7*sin(x3)-x4*0.1;
    x5dot = x7*cos(x3)-9.8-x5*0.1;
    x6dot = x8;
    x7dot = udf;
    x8dot = udtau;

    derivatives[ 0 ] = x1dot;
    derivatives[ 1 ] = x2dot;
    derivatives[ 2 ] = x3dot;
    derivatives[ 3 ] = x4dot;
    derivatives[ 4 ] = x5dot;
    derivatives[ 5 ] = x6dot;
    derivatives[ 6 ] = x7dot;
    derivatives[ 7 ] = x8dot;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    adouble x10 = initial_states[ 0 ];
    adouble x20 = initial_states[ 1 ];
    adouble x30 = initial_states[ 2 ];
    adouble x40 = initial_states[ 3 ];
    adouble x50 = initial_states[ 4 ];
    adouble x60 = initial_states[ 5 ];
    adouble x70 = initial_states[ 6 ];
    adouble x80 = initial_states[ 7 ];

    adouble x1f = final_states[ 0];
    adouble x2f = final_states[ 1];
    adouble x3f = final_states[ 2];
    adouble x4f = final_states[ 3];
    adouble x5f = final_states[ 4];
    adouble x6f = final_states[ 5];
    adouble x7f = final_states[ 6];
    adouble x8f = final_states[ 7];

    e[ 0 ] = x10;
    e[ 1 ] = x20;
    e[ 2 ] = x30;
    e[ 3 ] = x40;
    e[ 4 ] = x50;
    e[ 5 ] = x60;
    e[ 6 ] = x70;
    e[ 7 ] = x80;

    e[ 8 ] = x1f;
    e[ 9 ] = x2f;
    e[ 10 ] = x3f;
    e[ 11 ] = x4f;
    e[ 12 ] = x5f;
    e[ 13 ] = x6f;
    e[ 14 ] = x7f;
    e[ 15 ] = x8f;
}

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
    // No linkages as this is a single phase problem
}


void do_process()
{
    chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now();

    problem.phases(1).bounds.lower.events   	<< py0,  pz0,  phi0,  vy0,  vz0, omegax0, thrust0, tau0, pyf, pzf, phif, vyf,vzf,omegaxf,thrustf,tauf;
    problem.phases(1).bounds.upper.events   	<< py0,  pz0,  phi0,  vy0,  vz0, omegax0, thrust0, tau0, pyf, pzf, phif, vyf,vzf,omegaxf,thrustf,tauf;

    problem.phases(1).guess.controls = u;
    problem.phases(1).guess.states = x;
    problem.phases(1).guess.time = t;

    psopt(solution, problem, algorithm);
    chrono::time_point<chrono::steady_clock> end_time = chrono::steady_clock::now();
    int timeconsumption = chrono::duration_cast<chrono::milliseconds>(end_time - begin_time).count();
    cout << "time consume mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm: " << timeconsumption<< endl;

    if (solution.error_flag) return;

    x = solution.get_states_in_phase(1);
    u = solution.get_controls_in_phase(1);
    t = solution.get_time_in_phase(1);

}

void interpolation()
{
    int interpolationpoints=floor(t(0,pointnumber-1)*controlfreq)+1;
    int leftnode=0,rightnode=1;
    double leftweight=0,rightweight=0;
    for(int i=0;i<interpolationpoints;i++)
    {
        if(t(0,leftnode)<=(double)i/controlfreq && (double)i/controlfreq<=t(0,rightnode))
        {
            leftweight=(t(0,rightnode)-(double)i/controlfreq)/(t(0,rightnode)-t(0,leftnode));
            rightweight=((double)i/controlfreq-t(0,leftnode))/(t(0,rightnode)-t(0,leftnode));
            for(int j=0;j<8;j++)
            {
                x_refined(j,i)=x(j,leftnode)*leftweight+x(j,rightnode)*rightweight;
            }
        } else if((double)i/controlfreq>(t.row(0))(rightnode))
        {
            while((double)i/controlfreq>(t.row(0))(rightnode))
            {
                leftnode=rightnode;
                rightnode++;
            }
            leftweight=(t(0,rightnode)-(double)i/controlfreq)/(t(0,rightnode)-t(0,leftnode));
            rightweight=((double)i/controlfreq-t(0,leftnode))/(t(0,rightnode)-t(0,leftnode));
            for(int j=0;j<8;j++)
            {
                x_refined(j,i)=x(j,leftnode)*leftweight+x(j,rightnode)*rightweight;
            }
        }
    }
}

// rostopic subscribe
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
    ros::init(argc, argv, "minimizetime_control");
    ros::NodeHandle nh;
    ros::Rate rate(100);// it is the rate to check and process callback functions. the high rate is to avoid the latency to accept the pos and vel

    ros::Subscriber drone_imu_sub=nh.subscribe<sensor_msgs::Imu>("/mavros/imu/data",10,droneImu_cb);
    ros::Subscriber plane_rpy_sub = nh.subscribe<geometry_msgs::Vector3>("drone/current_rpy",10,dronerpy_cb);
    ros::Subscriber car_position_sub = nh.subscribe<nav_msgs::Odometry>("current_relative_postwist",10,relative_postwist_cb); //车的pos+twist
    ros::Subscriber plane_velocity_sub = nh.subscribe<geometry_msgs::TwistStamped>("mavros/local_position/velocity_local", 10, plane_vel_cb); //twist


    ros::Publisher controlstate_pub=nh.advertise<offb_posctl::controlstate>("ocp/control_state",10);
    ros::Publisher planeVel_pub=nh.advertise<geometry_msgs::Vector3>("filteredPlaneVel",10);

/////////////////////////////////// construct the ocp problem
    problem.name        		        = "minimize time Problem";
    problem.outfilename                 = "aggrecup.txt";
    problem.nphases   			        = 1;
    problem.nlinkages                   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   		= 8;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 16;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes                     << pointnumber;
//    problem.phases(1).nodes=(RowVectorXi(2)<<7,9).finished();
    psopt_level2_setup(problem, algorithm);


    problem.phases(1).bounds.lower.states   	<<  -1,  0.2,  -1.57, -10, -10,-10,  0,   -10;
    problem.phases(1).bounds.upper.states   	<<   3,  10,    1.57,  10,  10, 10,  19.6, 10;

    problem.phases(1).bounds.lower.controls 	<< -15, -15;
    problem.phases(1).bounds.upper.controls 	<<  15,  15;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 0.0;
    problem.phases(1).bounds.upper.EndTime      = 10.0;



    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae 		    = &dae;
    problem.events 		    = &events;
    problem.linkages		= &linkages;

    MatrixXd x0(8,pointnumber-1) ;
    MatrixXd u0(2,pointnumber-1) ;
    MatrixXd t0(1,pointnumber-1) ;

    x0.row(0) = linspace(0.0,1.0, pointnumber-1);
    x0.row(1) = linspace(0.0,1.0, pointnumber-1);
    x0.row(2) = linspace(0.0,1.0, pointnumber-1);
    x0.row(3) = linspace(0.0,1.0, pointnumber-1);
    x0.row(4) = linspace(0.0,1.0, pointnumber-1);
    x0.row(5) = linspace(0.0,1.0, pointnumber-1);
    x0.row(6) = linspace(0.0,1.0, pointnumber-1);
    x0.row(7) = linspace(0.0,1.0, pointnumber-1);

    u0.row(0) = linspace(0.0,1.0, pointnumber-1);
    u0.row(1) = linspace(0.0,1.0, pointnumber-1);

    t0.row(0) = linspace(0.0,2.0, pointnumber-1);

    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-3;
    algorithm.collocation_method          = "trapezoidal";
    algorithm.print_level                 =0;
/////////////////////////////////// construct the ocp problem

    problem.phases(1).bounds.lower.events   	<< py0,  pz0,  phi0,  vy0,  vz0, omegax0, thrust0, tau0, pyf, pzf, phif, vyf,vzf,omegaxf,thrustf,tauf;
    problem.phases(1).bounds.upper.events   	<< py0,  pz0,  phi0,  vy0,  vz0, omegax0, thrust0, tau0, pyf, pzf, phif, vyf,vzf,omegaxf,thrustf,tauf;

    problem.phases(1).guess.controls = u0;
    problem.phases(1).guess.states = x0;
    problem.phases(1).guess.time = t0;
    psopt(solution, problem, algorithm);

    if (solution.error_flag) return 0;

    x = solution.get_states_in_phase(1);
    u = solution.get_controls_in_phase(1);
    t = solution.get_time_in_phase(1);
    cout<<x<<endl;



    while (ros::ok())
    {
        ros::spinOnce();// to examine the queues of the callback functions once
        if(currentupdateflag)
        {
            do_process();
//            cout<<"solution.error_flag:  "<<solution.error_flag<<endl;
            if(!solution.error_flag)
            {
                interpolation();
                controlstate_msg.discrepointpersecond=controlfreq;
                controlstate_msg.inicounter=1;// the first value is current state, we should use at least next state
                controlstate_msg.arraylength=floor(t(0,pointnumber-1)*controlfreq)+1;// becasue the first point's index starts from 0.
                controlstate_msg.thrustarray.clear();
                controlstate_msg.phiarray.clear();
                controlstate_msg.tauarray.clear();
                controlstate_msg.stateYarray.clear();
                controlstate_msg.stateZarray.clear();
                controlstate_msg.stateVYarray.clear();
                controlstate_msg.stateVZarray.clear();
                for(int i=0;i<controlstate_msg.arraylength;i++)
                {
                    controlstate_msg.thrustarray.push_back(x_refined(6,i));
                    controlstate_msg.phiarray.push_back(x_refined(2,i));
                    controlstate_msg.thetaarray.push_back(0);
                    controlstate_msg.stateXarray.push_back(0);
                    controlstate_msg.stateYarray.push_back(x_refined(0,i));
                    controlstate_msg.stateZarray.push_back(x_refined(1,i));
                    controlstate_msg.stateVYarray.push_back(x_refined(3,i));
                    controlstate_msg.stateVZarray.push_back(x_refined(4,i));
                    controlstate_msg.tauarray.push_back(x_refined(7,i));
                }
//                cout<<"controlstate_msg.arraylength----!!!!:"<<controlstate_msg.arraylength<<"  time:"<<t(0,pointnumber-1)<<endl;
                controlstate_pub.publish(controlstate_msg);
//                cout<<"controlstate_msg.stateXarray[0]-----------fffffff:"<<controlstate_msg.stateYarray[0]<<endl;
            }
            currentupdateflag= false;
//            planeVel_pub.publish(filteredPlaneVelmsg);
        }
//        rate.sleep(); // comment this to achieve the maximize refresh frequency
    }

    return 0;
}
