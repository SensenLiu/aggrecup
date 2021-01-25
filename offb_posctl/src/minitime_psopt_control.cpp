//
// Created by lss on 2021/1/23.
//
#include "psopt.h"
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
    adouble udf = controls[ 0 ];
    adouble udtau = controls[ 1 ];
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



void do_process()
{
    double x10=0,x20=0.5,x30=0,x40=0,x50=0,x60=0,x70=9.8,x80=0;
    double x1f=2,x2f=2,x3f=1.57,x4f=0.2*sin(x3f),x5f=-0.2*cos(x3f),x6f=0.0,x7f=x70*cos(x3f),x8f=0;

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Brachistochrone Problem";

    problem.outfilename                 = "brac1.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 8;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 16;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes                     << 10;
//    problem.phases(1).nodes=(RowVectorXi(2)<<7,9).finished();
    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.phases(1).bounds.lower.states   	<<  -1,  0.5,  -1.57, -10, -10,-10,0,-10;
    problem.phases(1).bounds.upper.states   	<<   3, 10, 1.57,10,10,10,19.6,10;

    problem.phases(1).bounds.lower.controls 	<< -15,-15;
    problem.phases(1).bounds.upper.controls 	<< 15,15;

    problem.phases(1).bounds.lower.events   	<< x10,  x20,  x30,  x40,  x50,x60,x70,x80,x1f,x2f,x3f,x4f,x5f,x6f,x7f,x8f;
    problem.phases(1).bounds.upper.events   	<< x10,  x20,  x30,  x40,  x50,x60,x70,x80,x1f,x2f,x3f,x4f,x5f,x6f,x7f,x8f;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 0.0;
    problem.phases(1).bounds.upper.EndTime      = 10.0;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae 		= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;


////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////



    MatrixXd x0(8,9);
    MatrixXd u0(2,9);
    MatrixXd t0(1,9);

    x0.row(0) = linspace(0.0,1.0, 9);
    x0.row(1) = linspace(0.0,1.0, 9);
    x0.row(2) = linspace(0.0,1.0, 9);
    x0.row(3) = linspace(0.0,1.0, 9);
    x0.row(4) = linspace(0.0,1.0, 9);
    x0.row(5) = linspace(0.0,1.0, 9);
    x0.row(6) = linspace(0.0,1.0, 9);
    x0.row(7) = linspace(0.0,1.0, 9);

    u0.row(0) = linspace(0.0,1.0, 9);
    u0.row(1) = linspace(0.0,1.0, 9);

    t0.row(0) = linspace(0.0,2.0, 9);

    MatrixXd x 		= x0;
    MatrixXd u 		= u0;
    MatrixXd t 		= t0;


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-3;
//    algorithm.hessian			  = "exact";
    algorithm.collocation_method          = "trapezoidal";
//    algorithm.mesh_refinement             = "automatic";
    algorithm.print_level                 =0;
//    algorithm.ode_tolerance               =1.e-3;

    while(1)
    {

        x20=x20+0.01;

        problem.phases(1).bounds.lower.events   	<< x10,  x20,  x30,  x40,  x50,x60,x70,x80,x1f,x2f,x3f,x4f,x5f,x6f,x7f,x8f;
        problem.phases(1).bounds.upper.events   	<< x10,  x20,  x30,  x40,  x50,x60,x70,x80,x1f,x2f,x3f,x4f,x5f,x6f,x7f,x8f;

        problem.phases(1).guess.controls = u;
        problem.phases(1).guess.states = x;
        problem.phases(1).guess.time = t;
//        psopt_level2_setup(problem, algorithm);
////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////
        chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now();
        psopt(solution, problem, algorithm);
        chrono::time_point<chrono::steady_clock> end_time = chrono::steady_clock::now();
        int timeconsumption = chrono::duration_cast<chrono::milliseconds>(end_time - begin_time).count();
        cout << "time consume mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm: " << timeconsumption<< endl;
        if (solution.error_flag) exit(0);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

        x = solution.get_states_in_phase(1);
        u = solution.get_controls_in_phase(1);
        t = solution.get_time_in_phase(1);
//    MatrixXd H           = solution.get_dual_hamiltonian_in_phase(1);
//    MatrixXd lambda      = solution.get_dual_costates_in_phase(1);

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

//    Save(x,"x.dat");
//    Save(u,"u.dat");
//    Save(t,"t.dat");
//    Save(lambda,"p.dat");
        cout << x << endl;
        cout<< t<<endl;
    }

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
