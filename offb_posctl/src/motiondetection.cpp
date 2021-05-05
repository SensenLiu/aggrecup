//
// Created by sensenliu on 2021/4/29.
//
#include "ros/ros.h"
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include "ctime"
#include <math.h>
#include <list>


int main(int arc,char** arv)
{
    ros::init(arc,arv,"motiondetection");
    ros::NodeHandle nh;
    ros::Rate rate(30);
    ros::Publisher current_pos_pub=nh.advertise<nav_msgs::Odometry>("currenttarget_postwist",1);
    ros::Publisher targetactual_path_pub=nh.advertise<nav_msgs::Path>("currenttarget_path",1);

    nav_msgs::Odometry currenttarget_postwist_msg;
    nav_msgs::Path targetactual_path;
    geometry_msgs::PoseStamped targetactual_pose_stamped;
    double basevelocity_y=-0.5,basevelocity_z=0,currentvelocity_y=0,currentvelocity_z=0;
    double pos_y=4.0,pos_z=2.0;
    double noise_amplitued=0.05;
    int freq=30,counter=0;

    targetactual_path.header.stamp=ros::Time::now();
    targetactual_path.header.frame_id="ground_link";
    srand(time(0));
    while(ros::ok())
    {
        currentvelocity_y=basevelocity_y+noise_amplitued*(rand()%100+1-50)/100.0;
        currentvelocity_z=basevelocity_z+noise_amplitued*(rand()%100+1-50)/100.0;
//        std::cout<<"currentvelocity_y:"<<currentvelocity_y<<"  currentvelocity_z:"<<currentvelocity_z<<std::endl;
        pos_y=pos_y+1.0/freq*currentvelocity_y;
        pos_z=pos_z+1.0/freq*currentvelocity_z;
        currenttarget_postwist_msg.header.stamp=ros::Time::now();
        currenttarget_postwist_msg.pose.pose.position.x=0;
        currenttarget_postwist_msg.pose.pose.position.y=pos_y;
        currenttarget_postwist_msg.pose.pose.position.z=pos_z;
        currenttarget_postwist_msg.twist.twist.linear.x=0;
        currenttarget_postwist_msg.twist.twist.linear.y=currentvelocity_y;
        currenttarget_postwist_msg.twist.twist.linear.z=currentvelocity_z;

        targetactual_pose_stamped.pose.position.x=currenttarget_postwist_msg.pose.pose.position.x;
        targetactual_pose_stamped.pose.position.y=currenttarget_postwist_msg.pose.pose.position.y;
        targetactual_pose_stamped.pose.position.z=currenttarget_postwist_msg.pose.pose.position.z;
        targetactual_path.poses.push_back(targetactual_pose_stamped);
        if(targetactual_path.poses.size()>10*freq)
        {
            auto it=targetactual_path.poses.begin();
            targetactual_path.poses.erase(it);
        }

        current_pos_pub.publish(currenttarget_postwist_msg);
        targetactual_path_pub.publish(targetactual_path);
        counter++;

        ros::spinOnce();
        rate.sleep();
    }
    return 0;
}

