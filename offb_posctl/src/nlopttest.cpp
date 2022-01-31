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
#include "nlopt.hpp"

using namespace Eigen;
using namespace std;


typedef struct {
    double T, k, a0, u1, u2, x10, x20, x30, x40, x1tf, x2tf, x3tf, x4tf, lmd1, lmd2, lmd3, lmd4;
} my_constraint_data;

double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    my_constraint_data *d = reinterpret_cast<my_constraint_data*>(my_func_data);
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

//double myvconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data) {
//    my_constraint_data *d = (my_constraint_data *) data;
//    double a =1, b = 2;
//    if (!grad.empty()) {
//        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
//        grad[1] = -1.0;
//    }
//    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
//}


int main(int argc,char ** argv)
{
    ros::init(argc,argv,"motionprediction");
    ros::NodeHandle nh;
    ros::Rate rate(100);
    ros::Publisher predictionarray_pub=nh.advertise<offb_posctl::wallstate>("predictedwall_array",1, true);
    ros::Publisher predictionpath_pub=nh.advertise<nav_msgs::Path>("predictedwall_path",1);
    chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now();
    int timeconsumption = 0;
    int counte=0;
    while(ros::ok()&&counte<10) {
        counte++;
        int timenodenumber = 4;
        int constraintsnumber = 3;
        my_constraint_data data = {0.6, 1.5, 0.3, -1, 1, 8.7, 4.2, 0.96, 1.03, 11, 4.2, 1.02, 1.03, 1, 1, 1, 1};

        nlopt::opt opt(nlopt::LD_SLSQP, timenodenumber);
//    nlopt::opt opt(nlopt::LD_MMA, 2);
        std::vector<double> lb(timenodenumber, 0.0);
        std::vector<double> ub(timenodenumber, data.T);
        std::vector<double> constraints_tol(constraintsnumber, 1e-6);

        opt.set_lower_bounds(lb);
        opt.set_upper_bounds(ub);
        opt.set_min_objective(myvfunc, &data);
        opt.remove_inequality_constraints();
        opt.add_inequality_mconstraint(multiconstraint, NULL, constraints_tol);
        opt.set_xtol_rel(1e-4);
        opt.set_maxtime(0.01);

        std::vector<double> x(timenodenumber);
        x[0] = 0.0 / (timenodenumber + 1) * data.T;
        x[1] = 2.0 / (timenodenumber + 1) * data.T;
        x[2] = 3.0 / (timenodenumber + 1) * data.T;
        x[3] = 4.0 / (timenodenumber + 1) * data.T;
        cout << "initial guess: " << counte<<x[0] << "," << x[1] << "," << x[2] << "," << x[3] << endl;
        double minf;

        try {
            nlopt::result result = opt.optimize(x, minf);
            std::cout << "found minimum at f(" << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << ") = "
                      << std::setprecision(3) << minf << std::endl;
        }
        catch (std::exception &e) {
            std::cout << "nlopt failed: " << e.what() << std::endl;
        }
    }
    chrono::time_point<chrono::steady_clock> end_time = chrono::steady_clock::now();
    timeconsumption=chrono::duration_cast<chrono::milliseconds>(end_time - begin_time).count();
    cout <<"time acado consume~~...........~~~~~~~~mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm: "<<timeconsumption<< endl;

    while(ros::ok())
    {
        ros::spinOnce();
        rate.sleep();
    }
    return 0;
}