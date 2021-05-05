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
#include <nav_msgs/Path.h>

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
geometry_msgs::Point plane_expected_velocity; //车的零点和飞机的零点差3m，根据车的当前位置计算飞机位置
geometry_msgs::Point plane_expected_acceleration; //车的零点和飞机的零点差3m，根据车的当前位置计算飞机位置
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
offb_posctl::controlstate temp_controlstatearray_msg;

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
float az_ini = 0.0;
float ay_ini=0;
float phi_ini=0;
float omegax_ini=0;
float thrust_ini=9.8;
float tau_ini=0;
float t_end = 3;
double thrustforceacc = 0.0;
int pointnumber=150;// the number is almost always 20. It less, the accuracy won't be enough, if more, the time consumpiton will be too large.
int controlfreq=30;
int discretizedpointpersecond = (int)pointnumber/t_end;
int controlcounter=0;
int controlmode = 0;//0为最优控制,1为tractor
int quad_state=0;//0 climbhover, 1 AggressiveFly, 2 AdhesionPhase, 3 keep current state hover, 4 AdhesionSuccess
FILTER vy_Der(20),vz_Der(20),ay_Filter(20),az_Filter(20);


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
bool planstopflag=false;
bool startattitudecotrolflag=false;


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
    Quaterniond current_q(pose_drone_Imu.orientation.w, pose_drone_Imu.orientation.x, pose_drone_Imu.orientation.y, pose_drone_Imu.orientation.z);
    Quaterniond current_acc(0, pose_drone_Imu.linear_acceleration.x, pose_drone_Imu.linear_acceleration.y, pose_drone_Imu.linear_acceleration.z);
    Quaterniond accinword=current_q*current_acc*current_q.conjugate();

    ay_ini=accinword.y();
    az_ini=accinword.z()-9.8;
    ay_ini= ay_Filter.filter(ay_ini);
    az_ini= az_Filter.filter(az_ini);
}

void plane_vel_cb(const geometry_msgs::TwistStamped::ConstPtr &msg){
    vel_drone = *msg;
}

void car_pos_cb(const nav_msgs::Odometry::ConstPtr &msg) {
    pose_car_odom = *msg;

//    plane_expected_position.x = pose_car_odom.pose.pose.position.x; //-1 means axis difference
//    plane_expected_position.y = pose_car_odom.pose.pose.position.y; //-1 means axis difference
//    plane_expected_position.z = pose_car_odom.pose.pose.position.z + 0.5;
}
void plane_alt_cb(const std_msgs::Float64::ConstPtr &msg){
    plane_real_alt = *msg;
}
bool contstaterecieveflag= false;
int userfulpointcounter=1;
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
    ros::Subscriber controlstate_sub = nh.subscribe<offb_posctl::controlstate>("ocp/control_state", 1, controlstate_cb);
    // 【发布】飞机姿态/拉力信息 坐标系:NED系
    ros::Publisher ocplan_postwist_pub = nh.advertise<nav_msgs::Odometry>("ocplan_positiontwist", 1);//bvp计算的期望位置
    ros::Publisher ocplan_u_pub = nh.advertise<nav_msgs::Odometry>("ocplan_u", 1);//bvp计算的期望控制量(,没有转化的)
    ros::Publisher target_atti_thrust_pub = nh.advertise<mavros_msgs::AttitudeTarget>("mavros/setpoint_raw/attitude",
                                                                                      1);//发布给mavros的控制量(经过换算)
    ros::Publisher plane_rpy_pub = nh.advertise<geometry_msgs::Vector3>("drone/current_rpy", 1);//飞机当前的rpy角
    ros::Publisher current_relativepostwist_pub = nh.advertise<nav_msgs::Odometry>("current_relative_postwist",
                                                                                   1);//当前状态方程中的状态量,即相对量
    ros::Publisher targeterror_pub = nh.advertise<geometry_msgs::Vector3>("targeterror", 1);//目标误差
    ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("atucaltrajectory",1, true);

    // 频率 [30Hz]
    ros::Rate rate(controlfreq);   //50hz的频率发送/接收topic  ros与pixhawk之间,50Hz control frequency

    // log输出文件初始化
    logfile.open("/home/sensenliu/catkin_ws/src/gazebo_ros_learning/offb_posctl/data/pitch_log_hover1.csv", std::ios::out);
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

    Vector3d dronetowall_pos_inworld(0,0,0);
    Vector3d dronetowall_vel_inworld(0,0,0);
    Vector3d wall_zaxis(0,0,0);
    Vector3d wall_xaxis(1,0,0);
    Vector3d wall_yaxis(0,0,0);
    while (ros::ok()) {
        float cur_time = get_ros_time(begin_time_02);  // 当前时间
        if (planeupdateflag && pose_drone_odom.pose.pose.position.z>=0.2)   //订阅到飞机位置则flag为true，发布完相对位置的消息后flag置false
        {
            if(planstopflag==false)
            {
                px_ini = pose_drone_odom.pose.pose.position.x;
                py_ini = pose_drone_odom.pose.pose.position.y;
                pz_ini = pose_drone_odom.pose.pose.position.z;

                vx_ini = vel_drone.twist.linear.x;
                vy_ini = vel_drone.twist.linear.y;
                vz_ini = vel_drone.twist.linear.z;
//                ay_ini= ay_Filter.filter(vy_Der.derivation(vy_ini,cur_time));
//                az_ini= az_Filter.filter(vz_Der.derivation(vz_ini,cur_time));
//                ay_ini= ay_Filter.filter(ay_ini);
//                az_ini= az_Filter.filter(az_ini);
                ay_ini=min(max((double)ay_ini,-2*9.8),2*9.8);
                az_ini=min(max((double)az_ini,-9.8),9.8);
//                cout<<"ay_ini:"<<ay_ini<<" az_ini:"<<az_ini<<endl;


                pz_ini=min(max((double)pz_ini, 0.2),2.5);
                vz_ini=min(max((double)vz_ini, -5.0),5.0);
                vy_ini=min(max((double)vy_ini, -5.0),5.0);

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

//            std::cout << "py_ini:  " << py_ini << " pz_ini:  " << pz_ini << " phi_ini:  " << phi_ini << " vy_ini:  "
//                      << vy_ini <<" vz_ini: "<<vz_ini<<" omegax_ini: "<<omegax_ini<<" thrust_ini: "<<thrust_ini<<" tau_ini: "<<tau_ini<<std::endl;//输出,放到py文件中求解
            planeupdateflag = false;
        }

        ros::spinOnce();//刷新callback的消息
        cout<<"-------quad_state"<<quad_state<<endl;
        switch (quad_state) {
            case 0:
                plane_expected_position.z=plane_expected_position.z+ascentvel*cur_time;
                plane_expected_position.z=min(plane_expected_position.z,0.5);
                plane_expected_position.x=0;
                plane_expected_position.y=0;
                plane_expected_velocity.x=0;
                plane_expected_velocity.y=0;
                plane_expected_velocity.z=0;
                plane_expected_acceleration.x=0;
                plane_expected_acceleration.y=0;
                plane_expected_acceleration.z=0;

                pix_controller(cur_time);

                if(plane_expected_position.z>=0.5)
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
                cout<<"-------controlcounter"<<controlcounter<<endl;
                if(contstaterecieveflag)  //订阅到bvp计算的控制量则flag为true,用于起始时刻,还没算出bvp时
                {
                    lefnodeindex = controlcounter;
                    cout<<"-------lefnodeindex"<<lefnodeindex<<endl;
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
                    }
                    controlcounter++;
//                    if(controlcounter>=controlstatearray_msg.arraylength ||controlstatearray_msg.arraylength<=10)
//                    if((fabs(plane_expected_position.y-controlstatearray_msg.stateYarray[controlstatearray_msg.arraylength-1])<=0.1 && fabs(plane_expected_position.z-controlstatearray_msg.stateZarray[controlstatearray_msg.arraylength-1])<=0.1))
////                    if((controlstatearray_msg.arraylength/controlstatearray_msg.discrepointpersecond)<1.0)
//                    {
//                        planstopflag= true;
//                        ROS_ERROR_STREAM( "plane_expected_position.y: "<<plane_expected_position.y<<" plane_expected_position.z: "<<plane_expected_position.z<<" arraylength:"<<controlstatearray_msg.arraylength);
//                    }
                    if((controlstatearray_msg.arraylength-controlcounter)<0.15*controlstatearray_msg.discrepointpersecond)
                    {
                        startattitudecotrolflag=true;
                        ROS_ERROR_STREAM("startattitudecotrolflag:"<<startattitudecotrolflag<<" leftindex: "<<lefnodeindex<<" left_controlpoints: "<<controlstatearray_msg.arraylength-controlcounter<<" plane_expected_position.z: "<<plane_expected_position.z<<" pos_drone.pose.position.z: "<<pose_drone_odom.pose.pose.position.z);
//                        ROS_ERROR_STREAM( "(controlstatearray_msg.arraylength-controlcounter): "<<(controlstatearray_msg.arraylength-controlcounter)<<" lefttime:"<<((controlstatearray_msg.arraylength-controlcounter)/controlstatearray_msg.discrepointpersecond));
                    }
                    if(controlcounter>=controlstatearray_msg.arraylength)
                    {
                        ROS_ERROR_STREAM( "plane_expected_position.y:"<<plane_expected_position.y<<" plane_expected_position.z:"<<plane_expected_position.z<<" current z: "<<pose_drone_odom.pose.pose.position.z);
                        quad_state=2;
                        wall_zaxis << 0,controlstatearray_msg.stateAYarray[controlstatearray_msg.arraylength-1]/9.8,(controlstatearray_msg.stateAZarray[controlstatearray_msg.arraylength-1]+9.8)/9.8;
                        wall_yaxis=wall_zaxis.cross(wall_xaxis);

                    }
//                    planned_postwist_msg.pose.pose.orientation.x=-atan(controlstatearray_msg.stateAYarray[lefnodeindex]/(controlstatearray_msg.stateAZarray[lefnodeindex]+9.8));
                }
                pix_controller(cur_time);

                planned_postwist_msg.pose.pose.position.x=plane_expected_position.x;
                planned_postwist_msg.pose.pose.position.y=plane_expected_position.y;
                planned_postwist_msg.pose.pose.position.z=plane_expected_position.z;
                planned_postwist_msg.pose.pose.orientation.x=angle_target.x;
                planned_postwist_msg.twist.twist.angular.y=PIDVY.Output;
                planned_postwist_msg.twist.twist.angular.z=PIDVZ.Output;
                ROS_INFO_STREAM("Apporaching_thrust_target: "<< thrust_target<<" roll:"<<angle_target.x<<" pitch:"<<angle_target.y<<" yaw:"<<angle_target.z);
                break;
            case 2:
                startattitudecotrolflag=false;
                tempcounter++;
                if(tempcounter<=controlstatearray_msg.parabolictime*controlstatearray_msg.discrepointpersecond)
                {
                    tempgoalPx=pose_drone_odom.pose.pose.position.x;
                    tempgoalPy=pose_drone_odom.pose.pose.position.y;
                    tempgoalPz=pose_drone_odom.pose.pose.position.z;
                    cout<<"tempcounter: "<<tempcounter<<"  parabolictime*30: "<<controlstatearray_msg.parabolictime*30<<endl;
                    ROS_ERROR_STREAM("temp goal goal goal Px:"<<tempgoalPx<<" Py: "<<tempgoalPy<<" Pz: "<<tempgoalPz<<" pitch:"<<temp_angle.x);
                    thrust_target=param.THR_HOVER;
                    break;
                }
                tempCurrentPx=pose_drone_odom.pose.pose.position.x;
                tempCurrentPy=pose_drone_odom.pose.pose.position.y;
                tempCurrentPz=pose_drone_odom.pose.pose.position.z;
                dronetowall_pos_inworld << 0.0, tempCurrentPy-controlstatearray_msg.wall_y,tempCurrentPz-controlstatearray_msg.wall_z;
                dronetowall_vel_inworld << 0.0, vel_drone.twist.linear.y-controlstatearray_msg.wall_vy,vel_drone.twist.linear.z-controlstatearray_msg.wall_vz;
                ROS_ERROR_STREAM("current Px:"<<tempCurrentPx<<" Py: "<<tempCurrentPy<<" Pz: "<<tempCurrentPz<<" pitch:"<<temp_angle.x);
//                ROS_ERROR_STREAM("dronetowall_pos_inworld :"<<dronetowall_pos_inworld<<"wall_zaxis :"<<wall_zaxis<<" dronetowall_pos_inworld.dot(wall_zaxis) :"<<dronetowall_pos_inworld.dot(wall_zaxis)<<" dronetowall_vel_inworld.dot(wall_yaxis): "<<dronetowall_vel_inworld.dot(wall_yaxis));
                if(dronetowall_pos_inworld.dot(wall_zaxis)<-0.2||dronetowall_vel_inworld.dot(wall_yaxis)<-1.5)
//                if(tempCurrentPz<controlstatearray_msg.wall_z-0.25||tempCurrentPy<controlstatearray_msg.wall_y-0.20)
                {
                    quad_state=3;
                    tempcounter=0;
                    ROS_ERROR_STREAM("tempCurrentPz:"<<tempCurrentPz<<" controlstatearray_msg.wall_z: "<<controlstatearray_msg.wall_z<<" vel_read.z: "<<vel_drone.twist.linear.z);

                }
                if(tempcounter>=30&&(dronetowall_pos_inworld.dot(wall_zaxis)>=-0.2&&dronetowall_vel_inworld.dot(wall_yaxis)>=-1.5))
//                if(tempcounter>=30&&tempCurrentPz>=controlstatearray_msg.wall_z-0.15)
                {
                    quad_state=4;
                    tempcounter=0;
                    ROS_ERROR_STREAM("tempCurrentPz:"<<tempCurrentPz<<" controlstatearray_msg.wall_z: "<<controlstatearray_msg.wall_z<<" vel_read.z: "<<vel_drone.twist.linear.z);
                }
                orientation_target = euler2quaternion(-atan(controlstatearray_msg.stateAYarray[lefnodeindex]/(controlstatearray_msg.stateAZarray[lefnodeindex]+9.8)), 0, angle_target.z);
                thrust_target  = param.THR_HOVER*cos(-atan(controlstatearray_msg.stateAYarray[lefnodeindex]/(controlstatearray_msg.stateAZarray[lefnodeindex]+9.8)))*cos(0);   //目标推力值 to alleviate the gravity's component along the drone's z axis
                thrust_target  = max((double)thrust_target,param.THR_HOVER*0.5);   //目标推力值,只是用来保证提供扭矩，the drone is easy to fall freely and crash
//                thrust_target=impedancecontrol();// this may make the attitude out of control when the thrust is too low and cause the collapse in test without wall
                ROS_INFO_STREAM("Duringsuck_thrust_target: "<< thrust_target<<" roll:"<<-atan(controlstatearray_msg.stateAYarray[lefnodeindex]/(controlstatearray_msg.stateAZarray[lefnodeindex]+9.8))<<" pitch:"<<0<<" yaw:"<<angle_target.z);
                break;
            case 3:
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
                break;
        }
//        if(quad_state!=4&&quad_state!=2)
//        {
//            pix_controller(cur_time);
//        }
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

    float feedforwardcoefx=1/(1+fabs(acc_xd+acc_vxd)/(fabs(plane_expected_acceleration.x)+0.5));//0.5 is to avoid the zero of denominator
    float feedforwardcoefy=1/(1+fabs(acc_yd+acc_vyd)/(fabs(plane_expected_acceleration.y)+0.5));
    float feedforwardcoefz=1/(1+fabs(acc_zd+acc_vzd)/(fabs(plane_expected_acceleration.z)+0.5));
//    cout<<"feedforwardcoefx:"<<feedforwardcoefx<<" y:"<<feedforwardcoefy<<" z:"<<feedforwardcoefz<<endl;
    //计算输出
    if(startattitudecotrolflag==false)
    {
        PIDVX.Output=acc_xd+acc_vxd+feedforwardcoefx*plane_expected_acceleration.x;
        PIDVY.Output=acc_yd+acc_vyd+feedforwardcoefy*plane_expected_acceleration.y;
        PIDVZ.Output=acc_zd+acc_vzd+feedforwardcoefz*plane_expected_acceleration.z;
    }else{
        PIDVX.Output=plane_expected_acceleration.x;
        PIDVY.Output=plane_expected_acceleration.y;
        PIDVZ.Output=plane_expected_acceleration.z;
    }
//    cout<<"acc_p_error  planz_a:"<<plane_expected_acceleration.z<<"  vz_a:"<<acc_vzd<<"  z:"<<acc_zd<<" roll: "<<angle_receive.x<<" PIDVZ.Output: "<<PIDVZ.Output<<endl;

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
    logfile<<cur_time<<","<<rpy.y<<std::endl;

}