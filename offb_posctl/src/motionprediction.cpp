//
// Created by sensenliu on 2021/4/26.
//
#include "ros/ros.h"
#include <chrono>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/Vector3.h>
#include <geometry_msgs/Quaternion.h>
#include <geometry_msgs/Pose.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <mavros_msgs/State.h>
#include <mavros_msgs/AttitudeTarget.h>
#include <sensor_msgs/Imu.h>
#include <std_msgs/Bool.h>
#include <std_msgs/Float32.h>
#include "offb_posctl/wallstate.h"
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <thread>
#include <math.h>
#include <list>

using namespace Eigen;
using namespace std;

nav_msgs::Odometry current_wallstate;
nav_msgs::Path predictedwallpath_msg;
geometry_msgs::PoseStamped predictedwall_pose_stamped;

std::list<nav_msgs::Odometry> wallstate_list;
int wallstatelist_length=50;
offb_posctl::wallstate wallstate_msg;
double predict_freq=20.0;
double prediction_time=5.0;
int prediction_length=ceil(predict_freq*prediction_time);
Eigen::MatrixXf times_power_pos;
Eigen::MatrixXf times_power_vel;
bool wallstateupdateflag= false;
void wallpostwist_cb(const nav_msgs::Odometry::ConstPtr &msg)
{
    current_wallstate=*msg;
    wallstate_list.push_back(current_wallstate);
    if(wallstate_list.size()>wallstatelist_length)
    {
//        std::vector<nav_msgs::Odometry> ::iterator wallstate_iter = wallstate_list.begin();
//        wallstate_list.erase(wallstate_iter);
        wallstate_list.pop_front();
    }
    wallstateupdateflag= true;
}

void cal_time_matrix(Eigen::MatrixXf &times_power_pos,Eigen::MatrixXf &times_power_vel, int polyorder, int prediction_length)
{
    wallstate_msg.discrepointpersecond=predict_freq;
    wallstate_msg.arraylength=prediction_length;
    times_power_pos=Eigen::MatrixXf::Zero(prediction_length,polyorder+1);
    times_power_vel=Eigen::MatrixXf::Zero(prediction_length,polyorder+1);
    for(int i=0;i<times_power_pos.rows();i++)
    {
        for(int j=0;j<times_power_pos.cols();j++)
        {
            times_power_pos(i,j)=pow(i/predict_freq,polyorder-j);
            times_power_vel(i,j)=(polyorder-j)*pow(i/predict_freq,abs(polyorder-j-1));
//            cout<<"times_power_pos(i,j):"<<times_power_pos(i,j)<<" times_power_vel(i,j): "<<times_power_vel(i,j)<<endl;
        }
    }

}

bool motionprediction()
{
//    cout<<"times_power_pos.cols(): "<<times_power_pos.cols()<<endl;
    wallstateupdateflag= false;
    if(wallstate_list.size()!=wallstatelist_length)
    {
        return false;
    }
    double currenttime=ros::Time::now().toSec();
    Eigen::MatrixXf time_matrix(wallstatelist_length,times_power_pos.cols());
    Eigen::MatrixXf coeff_matrix(times_power_pos.cols(),3);// three dimensions pos,
    Eigen::MatrixXf output_matrix(wallstatelist_length,3);// three dimensions pos,
    Eigen::MatrixXf pos_predict_matrix(times_power_pos.rows(),3);
    Eigen::MatrixXf vel_predict_matrix(times_power_pos.rows(),3);

    time_matrix=Eigen::MatrixXf::Zero(wallstatelist_length,times_power_pos.cols());// must specify the dimensions 10 and 2 or the errors will occur
    coeff_matrix=Eigen::MatrixXf::Zero(times_power_pos.cols(),3);
    output_matrix=Eigen::MatrixXf::Zero(wallstatelist_length,3);
    pos_predict_matrix=Eigen::MatrixXf::Zero(times_power_pos.rows(),3);
    vel_predict_matrix=Eigen::MatrixXf::Zero(times_power_pos.rows(),3);
    auto wallstate_iterator=wallstate_list.begin();
    for(int i=0;i<wallstate_list.size();i++)//calculate the time matrix
    {
        for(int j=0;j<time_matrix.cols();j++)
        {
            time_matrix(i,j)=pow((*wallstate_iterator).header.stamp.toSec()-currenttime,times_power_pos.cols()-1-j);
//            cout<<"time_matrix(i,j):"<<time_matrix(i,j)<<endl;
        }
        output_matrix(i,0)=(*wallstate_iterator).pose.pose.position.x;
        output_matrix(i,1)=(*wallstate_iterator).pose.pose.position.y;
        output_matrix(i,2)=(*wallstate_iterator).pose.pose.position.z;
        wallstate_iterator++;
    }

    Eigen::JacobiSVD<Eigen::MatrixXf> svd(time_matrix);
    double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
    if (cond > 75)//max_cond_number: 自己设置的条件数阈值,我设置的是75
    {
        std::cout << "Matrix A is almost singular: "<<cond<<std::endl;
//        return false;
    }
    coeff_matrix=time_matrix.bdcSvd(ComputeThinU|ComputeThinV).solve(output_matrix);
//    cout<<"coeff_matrix: "<<coeff_matrix<<endl;
    wallstate_msg.header.stamp=ros::Time::now();
    wallstate_msg.arraylength=prediction_length;
    wallstate_msg.discrepointpersecond=(int)predict_freq;
    wallstate_msg.stateXarray.clear();
    wallstate_msg.stateYarray.clear();
    wallstate_msg.stateZarray.clear();
    wallstate_msg.stateVXarray.clear();
    wallstate_msg.stateVYarray.clear();
    wallstate_msg.stateVZarray.clear();

    predictedwallpath_msg.poses.clear();

    pos_predict_matrix=times_power_pos*coeff_matrix;
    vel_predict_matrix=times_power_vel*coeff_matrix;
//    cout<<"pos_predict_matrix: "<<pos_predict_matrix<<endl;
    for(int i=0;i<prediction_length;i++)
    {
//        cout<<"pos_predict_matrix(i,0): "<<pos_predict_matrix(i,0)<<endl;
        wallstate_msg.stateXarray.push_back(pos_predict_matrix(i,0));
        wallstate_msg.stateYarray.push_back(pos_predict_matrix(i,1));
        wallstate_msg.stateZarray.push_back(pos_predict_matrix(i,2));

        wallstate_msg.stateVXarray.push_back(vel_predict_matrix(i,0));
        wallstate_msg.stateVYarray.push_back(vel_predict_matrix(i,1));
        wallstate_msg.stateVZarray.push_back(vel_predict_matrix(i,2));
//        cout<<"wallstate_msg.stateVYarray:"<<wallstate_msg.stateVYarray[i]<<"wallstate_msg.stateVZarray:"<<wallstate_msg.stateVZarray[i]<<endl;


        predictedwall_pose_stamped.pose.position.x=wallstate_msg.stateXarray[i];
        predictedwall_pose_stamped.pose.position.y=wallstate_msg.stateYarray[i];
        predictedwall_pose_stamped.pose.position.z=wallstate_msg.stateZarray[i];
        predictedwallpath_msg.poses.push_back(predictedwall_pose_stamped);
    }
    return true;

}
int main(int argc,char ** argv)
{
    ros::init(argc,argv,"motionprediction");
    ros::NodeHandle nh;
    ros::Rate rate(100);
    ros::Subscriber car_position_sub = nh.subscribe<nav_msgs::Odometry>("currenttarget_postwist",10,wallpostwist_cb); //车的pos+twist
    ros::Publisher predictionarray_pub=nh.advertise<offb_posctl::wallstate>("predictedwall_array",1, true);
    ros::Publisher predictionpath_pub=nh.advertise<nav_msgs::Path>("predictedwall_path",1);

    predictedwallpath_msg.header.stamp=ros::Time::now();
    predictedwallpath_msg.header.frame_id="ground_link";

    cal_time_matrix(times_power_pos,times_power_vel, 1, prediction_length);
    while(ros::ok())
    {
        if(wallstateupdateflag)
        {
            if(motionprediction())
            {
                predictionarray_pub.publish(wallstate_msg);
                predictionpath_pub.publish(predictedwallpath_msg);
            }
        }
        ros::spinOnce();
        rate.sleep();
    }
    return 0;
}