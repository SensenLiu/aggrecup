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


#include "ros/ros.h"
#include "std_msgs/Float32.h"
#include <chrono>
#include "iomanip"
#include <thread>
#include <unistd.h>
#include <mutex>
#include "offb_posctl/controlstate.h"

using namespace Eigen;//释放eigen命名空间 矩阵库
using namespace std;

///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>全 局 变 量<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

mavros_msgs::State current_state;           //无人机当前状态
nav_msgs::Odometry pose_drone_odom;       //读入的无人机drone当前位置，x，y，z+姿态
sensor_msgs::Imu pose_drone_Imu;       //读入的无人机drone当前位置，x，y，z+姿态
nav_msgs::Odometry pose_car_odom;       //读入car当前位置
geometry_msgs::TwistStamped vel_drone;      //读入的无人机当前速度 线速度+角速度
geometry_msgs::Quaternion orientation_target;   //发给无人机的姿态指令  四元数
geometry_msgs::Vector3 angle_target;   //欧拉角
geometry_msgs::Vector3 vel_target;   //期望速度
geometry_msgs::Point plane_expected_position; //车的零点和飞机的零点差3m，根据车的当前位置计算飞机位置
geometry_msgs::PoseStamped target_attitude;  //1
mavros_msgs::Thrust target_thrust_msg; //1循环
std_msgs::Float64 plane_real_alt; //control前
mavros_msgs::AttitudeTarget target_atti_thrust_msg; //最终发布的消息 油门+角度
mavros_msgs::AttitudeTarget base_atti_thrust_msg; //最终发布的消息 油门+角度
//mavros_msgs::PositionTarget target_pos_msg;
nav_msgs::Odometry  planned_postwist_msg;
nav_msgs::Odometry planned_u_msg;
nav_msgs::Odometry current_relativepostwist_msg;
geometry_msgs::Vector3 targeterror_msg;
geometry_msgs::Vector3 temp_angle;
geometry_msgs::Vector3 rpy;
offb_posctl::controlstate controlstatearray_msg;

float thrust_target;        //期望推力
float Yaw_Init;
float Yaw_Locked = 0;           //锁定的偏航角(一般锁定为0)
bool got_initial_point = false;
PID PIDVX, PIDVY, PIDVZ;    //声明PID类
Parameter param;
std::ofstream logfile;

///for psopt
float px_ini = -3.0;
float pz_ini = 0;
float py_ini=0;
float vx_ini = -0.1;
float vz_ini = 0.0;
float vy_ini=0;
float phi_ini=0;
float omegax_ini=0;
float thrust_ini=9.8;
float tau_ini=0;
FILTER derivation_omegax(3);
float t_end = 3;
double thrustforceacc = 0.0;
int pointnumber=150;// the number is almost always 20. It less, the accuracy won't be enough, if more, the time consumpiton will be too large.
int controlfreq=50;
int discretizedpointpersecond = (int)pointnumber/t_end;
int controlcounter=0;
//int controltime = 15; //取的控制序列点数
int controlmode = 0;//0为最优控制,1为tractor
bool adjust_flag = true;
double euler_anlge_limit=0.314;


///for pitch compensation
float comp_integrate,comp_last_error;
float comp_kp = 0,comp_ki = 0,comp_kd = 0;


/// for tractor
std::mutex mtx;
double amp=0.5;
float sinrate=3;
double ocpPitch=0;
double ocpRoll=0;
double ax=0.0,az=0.0,ay=0.0;
double axp=0.0,azp=0.0,ayp=0.0;
double axv=0.0,azv=0.0,ayv=0.0;

///virtual car pos_twist
float virtual_pos_x = 0;
float virtual_vel_x = 0;

///
double bvp_restrict_ax, bvp_restrict_az;// bvp acc_x&z restriction
double bvp_feedback_gain_x[2] = {0.5,0.3}, bvp_feedback_gain_z[2] = {0.5,0.3}; ///bvp pos&vel feedback gain
double bvp_acc_x, bvp_acc_y,bvp_acc_z, bvp_ocpPitch, bvp_thrustforceacc; ///bvp feedback acc desired
///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>声 明 函 数<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

//欧拉角转四元数
geometry_msgs::Quaternion euler2quaternion(float roll, float pitch, float yaw);//geometry_msgs的Quaternion类型的函数
geometry_msgs::Vector3 quaternion2euler(float x, float y, float z, float w);

float get_ros_time(ros::Time time_begin);                                            //获取ros当前时间
int pix_controller(float cur_time);
void vector3dLimit(Vector3d &v, double limit) ; ///limit should be positive
Vector3d vectorElementMultiply(Vector3d v1, Vector3d v2);
void tractor_controller(float time);
//int pix_controller(int cur_time);
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
}
void plane_imu_cb(const sensor_msgs::Imu::ConstPtr &msg){
    pose_drone_Imu = *msg;
}

void plane_vel_cb(const geometry_msgs::TwistStamped::ConstPtr &msg){
    vel_drone = *msg;
}

void car_pos_cb(const nav_msgs::Odometry::ConstPtr &msg) {
    pose_car_odom = *msg;

    plane_expected_position.x = pose_car_odom.pose.pose.position.x; //-1 means axis difference
    plane_expected_position.y = pose_car_odom.pose.pose.position.y; //-1 means axis difference
    plane_expected_position.z = pose_car_odom.pose.pose.position.z + 0.5;
}
void plane_alt_cb(const std_msgs::Float64::ConstPtr &msg){
    plane_real_alt = *msg;
}
bool contstaterecieveflag= false;
void controlstate_cb(const offb_posctl::controlstate::ConstPtr &msg)
{
    controlstatearray_msg = *msg;
    controlcounter = controlstatearray_msg.inicounter;
    if(contstaterecieveflag == false)//第一次回调时初始化，之后这个flag一直是true
    {
        contstaterecieveflag= true;
        discretizedpointpersecond=controlstatearray_msg.discrepointpersecond;
    }

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
    ros::Subscriber plane_alt_sub = nh.subscribe<std_msgs::Float64>("mavros/global_position/rel_alt", 1, plane_alt_cb);
    ros::Subscriber controlstate_sub = nh.subscribe<offb_posctl::controlstate>("bvp_controlstate", 1, controlstate_cb);
    // 【发布】飞机姿态/拉力信息 坐标系:NED系
    ros::Publisher ocplan_postwist_pub = nh.advertise<nav_msgs::Odometry>("ocplan_positiontwist", 1);//bvp计算的期望位置
    ros::Publisher ocplan_u_pub = nh.advertise<nav_msgs::Odometry>("ocplan_u", 1);//bvp计算的期望控制量(,没有转化的)
    ros::Publisher target_atti_thrust_pub = nh.advertise<mavros_msgs::AttitudeTarget>("mavros/setpoint_raw/attitude",
                                                                                      1);//发布给mavros的控制量(经过换算)
    ros::Publisher plane_rpy_pub = nh.advertise<geometry_msgs::Vector3>("drone/current_rpy", 1);//飞机当前的rpy角
    ros::Publisher current_relativepostwist_pub = nh.advertise<nav_msgs::Odometry>("current_relative_postwist",
                                                                                   1);//当前状态方程中的状态量,即相对量
    ros::Publisher targeterror_pub = nh.advertise<geometry_msgs::Vector3>("targeterror", 1);//目标误差
    // 频率 [30Hz]
    ros::Rate rate(controlfreq);   //50hz的频率发送/接收topic  ros与pixhawk之间,50Hz control frequency

    // log输出文件初始化
    logfile.open("/home/sensenliu/catkin_ws/src/gazebo_ros_learning/offb_posctl/log/pitch_log_hover1.csv", std::ios::out);
    if (!logfile.is_open()) {
        ROS_ERROR("log to file error!");
//        return 0;
    }

    // 读取PID参数
    std::string paraadr("/home/sensenliu/catkin_ws/src/gazebo_ros_learning/offb_posctl/src/param");
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
    target_atti_thrust_msg.thrust = 0.65; //65%的油门 50%与重力平衡即悬停

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
                }
                last_request = ros::Time::now();
            }
        }
        target_atti_thrust_pub.publish(target_atti_thrust_msg);
        if (plane_real_alt.data > 0.7) {
            ROS_INFO("plane takeoff !");
            break;
        }
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
    int quad_state=0;//0 climbhover, 1 AggressiveFly, 2 AdhesionPhase, 3 keep current state hover, 4 AdhesionSuccess
    int tempcounter=0;
    float ascentvel=15;
    float tempCurrentPx=0,tempCurrentPy=0,tempCurrentPz=0;
    int lefnodeindex = 0;
    int rightnodeindex = 0;

///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>主  循  环<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    while (ros::ok()) {

        ros::spinOnce();//刷新callback的消息

        float cur_time = get_ros_time(begin_time_02);  // 当前时间
        switch (quad_state) {
            case 0:
                plane_expected_position.z=plane_expected_position.z+ascentvel*cur_time;
                plane_expected_position.z=min(plane_expected_position.z,0.5);
                plane_expected_position.x=0;
                plane_expected_position.y=0;
                pix_controller(cur_time);

                if(plane_expected_position.z>=0.5)
                {
                    tempcounter++;
                    if(tempcounter>=150)
                    {
                        quad_state=1;
                        controlcounter=1;
                        tempcounter=0;
                    }
                }
                break;
            case 1:
                if(contstaterecieveflag)  //订阅到bvp计算的控制量则flag为true,用于起始时刻,还没算出bvp时
                {
                    lefnodeindex = controlcounter;
                    if (lefnodeindex+1 < controlstatearray_msg.arraylength)
                    {
                        ///compensation
                        ocpRoll = controlstatearray_msg.phiarray[lefnodeindex];
                        thrustforceacc = controlstatearray_msg.thrustarray[lefnodeindex];

                        plane_expected_position.z=controlstatearray_msg.stateZarray[lefnodeindex];
                        plane_expected_position.x=controlstatearray_msg.stateXarray[lefnodeindex];
                        plane_expected_position.y=controlstatearray_msg.stateYarray[lefnodeindex];
                        pix_controller(cur_time);

                        if (lefnodeindex <= 6)
                        {
                            orientation_target = euler2quaternion(ocpRoll, angle_target.y, angle_target.z);
                            thrust_target = (float) (param.hoverthrust) * thrustforceacc / 9.8;
                        }
                    }
                    controlcounter++;
                    if(controlcounter>=controlstatearray_msg.arraylength ||controlstatearray_msg.arraylength<=10)
                    {
                        quad_state=2;
                    }
                } else{
                    pix_controller(cur_time);
                }
                break;
            case 2:
                tempCurrentPx=pose_drone_odom.pose.pose.position.x;
                tempCurrentPy=pose_drone_odom.pose.pose.position.y;
                tempCurrentPz=pose_drone_odom.pose.pose.position.z;;
                tempcounter++;
                if(tempCurrentPz<controlstatearray_msg.stateZarray[controlstatearray_msg.arraylength-1]-0.15)
                {
                    quad_state=3;
                    tempcounter=0;
                }
                if(tempcounter>=100&&tempCurrentPz>=controlstatearray_msg.stateZarray[controlstatearray_msg.arraylength-1]-0.15)
                {
                    quad_state=4;
                    tempcounter=0;
                }
                orientation_target = euler2quaternion(controlstatearray_msg.phiarray[controlstatearray_msg.arraylength-1], controlstatearray_msg.thetaarray[controlstatearray_msg.arraylength-1], angle_target.z);
                thrust_target  = param.hoverthrust*cos(controlstatearray_msg.phiarray[controlstatearray_msg.arraylength-1])*cos(controlstatearray_msg.thetaarray[controlstatearray_msg.arraylength-1]);   //目标推力值 to alleviate the gravity's component along the drone's z axis
//                thrust_target  = 0.2;   //目标推力值,只是用来保证提供扭矩，the drone is easy to fall freely and crash
                ROS_INFO_STREAM("Duringsuck_thrust_target: "<< thrust_target<<" roll:"<<controlstatearray_msg.phiarray[controlstatearray_msg.arraylength-1]<<" pitch:"<<controlstatearray_msg.thetaarray[controlstatearray_msg.arraylength-1]<<" yaw:"<<angle_target.z);
                break;
                break;
            case 3:
                plane_expected_position.x=tempCurrentPx;
                plane_expected_position.y=tempCurrentPy;
                plane_expected_position.z=tempCurrentPz;
                pix_controller(cur_time);
                break;
            case 4:
                orientation_target = euler2quaternion(controlstatearray_msg.phiarray[controlstatearray_msg.arraylength-1], controlstatearray_msg.thetaarray[controlstatearray_msg.arraylength-1], angle_target.z);
                thrust_target  = 0;   //目标推力值
                ROS_INFO_STREAM("Sucksuccess_thrust_target: "<< thrust_target<<" roll:"<<angle_target.x<<" pitch:"<<angle_target.y<<" yaw:"<<angle_target.z);
                break;
            default:
                break;
        }
//        if(quad_state!=4&&quad_state!=2)
//        {
//            pix_controller(cur_time);
//        }
        cout<<"refpos x y z: "<<plane_expected_position.x<<"   y:"<<plane_expected_position.y<<"  z:"<<plane_expected_position.z<<endl;

        if (planeupdateflag && pose_drone_odom.pose.pose.position.z>=0.2)   //订阅到飞机位置则flag为true，发布完相对位置的消息后flag置false
        {
            omegax_ini=min(max(pose_drone_odom.twist.twist.angular.x,-10.0),10.0);

            tau_ini=derivation_omegax.derivation(cur_time,pose_drone_odom.twist.twist.angular.x);
            tau_ini=min(max((double)tau_ini,-10.0),10.0);

            phi_ini=temp_angle.x;

            px_ini = pose_drone_odom.pose.pose.position.x;
            pz_ini = max(pose_drone_odom.pose.pose.position.z,0.2);
            py_ini = pose_drone_odom.pose.pose.position.y;
            vx_ini = vel_drone.twist.linear.x;
            vz_ini = vel_drone.twist.linear.z;
            vy_ini = vel_drone.twist.linear.y;

            thrust_ini=thrust_target/param.hoverthrust*9.8;
            thrust_ini=min((double)thrust_ini,19.6);


            std::cout << "px_ini:  " << px_ini << "pz_ini:  " << pz_ini << "vx_ini:  " << vx_ini << "vz_ini:  "
                      << vz_ini << std::endl;//输出,放到py文件中求解
            cout<<"va_ini:"<<vel_drone.twist.linear.x<<endl;
            current_relativepostwist_msg.pose.pose.position.x = px_ini;
            current_relativepostwist_msg.pose.pose.position.z = pz_ini;
            current_relativepostwist_msg.pose.pose.position.y = py_ini;
            current_relativepostwist_msg.twist.twist.linear.x = vx_ini;
            current_relativepostwist_msg.twist.twist.linear.y = vy_ini;
            current_relativepostwist_msg.twist.twist.linear.z = vz_ini;
            current_relativepostwist_msg.pose.pose.orientation.x=phi_ini;
            current_relativepostwist_msg.pose.pose.orientation.y=omegax_ini;
            current_relativepostwist_msg.pose.pose.orientation.z=thrust_ini;
            current_relativepostwist_msg.pose.pose.orientation.w=tau_ini;
            current_relativepostwist_msg.header.stamp = pose_drone_odom.header.stamp;
            current_relativepostwist_pub.publish(current_relativepostwist_msg);
            planeupdateflag = false;
        }

            ///publish plane current rpy
            temp_angle = quaternion2euler(pose_drone_odom.pose.pose.orientation.x,
                                          pose_drone_odom.pose.pose.orientation.y,
                                          pose_drone_odom.pose.pose.orientation.z,
                                          pose_drone_odom.pose.pose.orientation.w);//欧拉角
            rpy.y = temp_angle.y;
            rpy.x = temp_angle.x;
            rpy.z = temp_angle.z;
            plane_rpy_pub.publish(rpy);


            ///publish thrust & orientation
            std::cout << "thrust_target: " << thrust_target << std::endl;
            target_atti_thrust_msg.header.stamp = ros::Time::now();
            target_atti_thrust_msg.orientation = orientation_target;
            target_atti_thrust_msg.thrust = thrust_target;
            target_atti_thrust_pub.publish(target_atti_thrust_msg);

            ///publish planned pos&twist in x&z
            planned_postwist_msg.header.stamp = ros::Time::now();
            ocplan_postwist_pub.publish(planned_postwist_msg);

            ///publish planned thrust & pitch
            ocplan_u_pub.publish(planned_u_msg);

            ///publish targeterror_msg
            targeterror_msg.x = px_ini + 3;
            targeterror_msg.z = pz_ini + px_ini * rpy.y;
            targeterror_msg.y = py_ini;//pub time consumption
            targeterror_pub.publish(targeterror_msg);
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

void tractor_controller(float time)
{
    static Vector3d z_w_norm(0, 0, 1.0);
    Vector3d p_error;
    Vector3d v_error;
    Vector3d planned_a;
    switch (controlmode)
    {
        case 0 :{
            p_error << (planned_postwist_msg.pose.pose.position.x - px_ini), (pose_car_odom.pose.pose.position.y-py_ini), (planned_postwist_msg.pose.pose.position.z-pz_ini);
            v_error << (planned_postwist_msg.twist.twist.linear.x - vx_ini), (pose_car_odom.twist.twist.linear.y-vy_ini), (planned_postwist_msg.twist.twist.linear.z-vz_ini);
            planned_a << sin(ocpPitch)*thrustforceacc, 0, thrustforceacc-9.8;
            break;
        }
        case 1 :{
            p_error << (planned_postwist_msg.pose.pose.position.x - px_ini), (pose_car_odom.pose.pose.position.y-py_ini), (planned_postwist_msg.pose.pose.position.z-pz_ini);
            v_error << (planned_postwist_msg.twist.twist.linear.x - vx_ini), (pose_car_odom.twist.twist.linear.y-vy_ini), (planned_postwist_msg.twist.twist.linear.z-vz_ini);
            planned_a << -amp*sinrate*sinrate*sin(sinrate*time), 0, 0;
            break;
        }
    }

    static Vector3d p_error_last;
    static Vector3d v_error_last;
    static Vector3d p_error_accumulate;
    static Vector3d v_error_accumulate;
    static bool if_init = true;
    static Vector3d position_error_p(param.txp_p,param.typ_p,param.tzp_p);
    static Vector3d position_error_d(param.txp_d,param.typ_d,param.tzp_d);
    static Vector3d position_error_i(param.txp_i,param.typ_i,param.tzp_i);
    static Vector3d velocity_error_p(param.txv_p,param.tyv_p,param.tzv_p);
    static Vector3d velocity_error_d(param.txv_d,param.tyv_d,param.tzv_d);
    static Vector3d velocity_error_i(param.txv_i,param.tyv_i,param.tzv_i);

    if(if_init){
        if_init = false;
        p_error_last = p_error;
        v_error_last = v_error;
        p_error_accumulate = p_error;
        v_error_accumulate = v_error;
        return;
    }

    /**Core code**/
    Vector3d delt_p_error = p_error - p_error_last;
    Vector3d delt_v_error = v_error - v_error_last;

    p_error_accumulate += p_error;
    v_error_accumulate += v_error;
    vector3dLimit(p_error_accumulate, 0.6);
    vector3dLimit(v_error_accumulate, 0.5);

    Vector3d a_fb =   /// PID
            vectorElementMultiply(p_error, position_error_p) + vectorElementMultiply(v_error, velocity_error_p) +
            vectorElementMultiply(delt_p_error, position_error_d) + vectorElementMultiply(delt_v_error, velocity_error_d) +
            vectorElementMultiply(p_error_accumulate, position_error_i) + vectorElementMultiply(v_error_accumulate, velocity_error_i);

    ax=a_fb(0);
    az=a_fb(2);
    ay=a_fb(1);
    axp=(vectorElementMultiply(p_error, position_error_p)+vectorElementMultiply(delt_p_error, position_error_d)+vectorElementMultiply(p_error_accumulate, position_error_i))(0);
    azp=(vectorElementMultiply(p_error, position_error_p)+vectorElementMultiply(delt_p_error, position_error_d)+vectorElementMultiply(p_error_accumulate, position_error_i))(2);
    ayp=(vectorElementMultiply(p_error, position_error_p)+vectorElementMultiply(delt_p_error, position_error_d)+vectorElementMultiply(p_error_accumulate, position_error_i))(1);
    axv=(vectorElementMultiply(v_error, velocity_error_p)+vectorElementMultiply(delt_v_error, velocity_error_d)+ vectorElementMultiply(v_error_accumulate, velocity_error_i))(0);
    azv=(vectorElementMultiply(v_error, velocity_error_p)+vectorElementMultiply(delt_v_error, velocity_error_d)+ vectorElementMultiply(v_error_accumulate, velocity_error_i))(2);
    ayv=(vectorElementMultiply(v_error, velocity_error_p)+vectorElementMultiply(delt_v_error, velocity_error_d)+ vectorElementMultiply(v_error_accumulate, velocity_error_i))(1);

    p_error_last = p_error;
    v_error_last = v_error;


    Vector3d a_des = a_fb + planned_a + 9.8 * z_w_norm;
    Vector3d att_des_norm = a_des / a_des.norm();
    ///quaternion way to determine attitude
//    Quaterniond att_des_q = Quaterniond::FromTwoVectors(z_w_norm, att_des_norm);
//    //add yaw
//    Quaterniond yaw_quat(cos(0/2.0), att_des_norm(0)*sin(0/2.0),
//                         att_des_norm(1)*sin(0/2.0),att_des_norm(2)*sin(0/2.0));
//    att_des_q = yaw_quat * att_des_q;
//
//    //Calculate thrust
//    orientation_target.x = att_des_q.x();
//    orientation_target.y = att_des_q.y();
//    orientation_target.z = att_des_q.z();
//    orientation_target.w = att_des_q.w();
    ///quaternion way to determine attitude

    angle_target.x = asin(-a_des(1)/a_des.norm());
    angle_target.y = atan(a_des(0)/a_des(2));
    angle_target.z = Yaw_Init;
    orientation_target = euler2quaternion(angle_target.x, angle_target.y, angle_target.z);

//    thrust_target  = (float)a_des.norm() /9.8*(base_atti_thrust_msg.thrust);   //目标推力值

    temp_angle = quaternion2euler(pose_drone_odom.pose.pose.orientation.x, pose_drone_odom.pose.pose.orientation.y,
                                  pose_drone_odom.pose.pose.orientation.z, \
        pose_drone_odom.pose.pose.orientation.w);
    thrust_target  = (float)(a_des(0)*sin(temp_angle.y)*cos(temp_angle.x)-a_des(1)*sin(temp_angle.x)+a_des(2)*cos(temp_angle.y)*cos(temp_angle.x)) /9.8*(base_atti_thrust_msg.thrust);   //目标推力值
}

int pix_controller(float cur_time)
{
//位 置 环
    //计算误差
    float error_x = plane_expected_position.x - pose_drone_odom.pose.pose.position.x;
    float error_y = plane_expected_position.y - pose_drone_odom.pose.pose.position.y;
    float error_z = plane_expected_position.z - plane_real_alt.data;
//    std::cout << "error: x：" << error_x << "\ty：" << error_y << "\tz：" << error_z << std::endl;
    //计算指定速度误差
    float vel_xd = param.x_p * error_x;
    float vel_yd = param.y_p * error_y;
    float vel_zd = param.z_p * error_z;
    vel_target.x = vel_xd;
    vel_target.y = vel_yd;
    vel_target.z = vel_zd;

//速 度 环
    //积分标志位.未进入OFFBOARD时,不累积积分项;进入OFFBOARD时,开始积分.
    PIDVX.start_intergrate_flag = true;
    PIDVY.start_intergrate_flag = true;
    PIDVZ.start_intergrate_flag = true;
    if(got_initial_point == false){
        PIDVX.start_intergrate_flag = false;
        PIDVY.start_intergrate_flag = false;
        PIDVZ.start_intergrate_flag = false;
    }
    //计算误差
    float error_vx = vel_xd - vel_drone.twist.linear.x;
    float error_vy = vel_yd - vel_drone.twist.linear.y;
    float error_vz = vel_zd - vel_drone.twist.linear.z;
    //传递误差
    PIDVX.add_error(error_vx, cur_time); //把error放到list中
    PIDVY.add_error(error_vy, cur_time);
    PIDVZ.add_error(error_vz, cur_time);
    //计算输出
    PIDVX.pid_output();
    PIDVY.pid_output();
    PIDVZ.pid_output();

//    Matrix2f A_yaw;
//    A_yaw << sin(Yaw_Locked), cos(Yaw_Locked),
//            -cos(Yaw_Locked), sin(Yaw_Locked);
//    Vector2f mat_temp(PIDVX.Output,PIDVY.Output);       //赋值到期望推力和姿态 x是前后，y是左右
//    Vector2f euler_temp= 1/9.8 * A_yaw.inverse() * mat_temp;
//    angle_target.x = euler_temp[0];
//    angle_target.y = euler_temp[1];
//    std::cout << " PIDVX.pid_output(): " << PIDVX.Output << "\tangle_target.y: " << angle_target.y << std::endl;
////    angle_target.z = Yaw_Locked + Yaw_Init;
//    angle_target.z = Yaw_Init;

    angle_target.x = asin(-PIDVY.Output/sqrt(pow(PIDVX.Output,2)+pow(PIDVY.Output,2)+pow(PIDVZ.Output+9.8,2)));
    angle_target.y = atan(PIDVX.Output/(PIDVZ.Output+9.8));
    angle_target.z = Yaw_Init;

    orientation_target = euler2quaternion(angle_target.x, angle_target.y, angle_target.z);
//    thrust_target = (float)(0.05 * (9.8 + PIDVZ.Output));   //目标推力值
    thrust_target  = (float)sqrt(pow(PIDVX.Output,2)+pow(PIDVY.Output,2)+pow(PIDVZ.Output+9.8,2))/9.8*(param.hoverthrust);   //目标推力值

//    std::cout << "PIDVZ.OUTPUT:  " << PIDVZ.Output << std::endl;
//    std::cout << "thrust_target:  " << thrust_target << std::endl;

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
    logfile<<cur_time<<","<<rpy.y<<std::endl;

}