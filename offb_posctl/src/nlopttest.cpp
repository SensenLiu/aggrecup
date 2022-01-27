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
    double a, b;
} my_constraint_data;

my_constraint_data data[2] = { {2,0}, {-1,1} };

double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty()) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
}


double myvconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    my_constraint_data *d = reinterpret_cast<my_constraint_data*>(data);
    double a = d->a, b = d->b;
    if (!grad.empty()) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}

void multiconstraint(unsigned m, double *result, unsigned n, const double *x, double *gradient, void *func_data) // ref:https://blog.csdn.net/weixin_33937499/article/details/93488349 and https://nlopt.readthedocs.io/en/latest/NLopt_Reference/
{
    my_constraint_data *d = reinterpret_cast<my_constraint_data*>(func_data);
    double a1 = d[0].a, b1 = d[0].b;
    double a2 = d[1].a, b2 = d[1].b;
    if(gradient){
//        gradient[0]= 1.0;
//        gradient[1]= -1.0;
//        gradient[2]= 0.0;
//        gradient[3]= 0.0;
//
//        gradient[4]= 0.0;
//        gradient[5]= 1.0;
//        gradient[6]= -1.0;
//        gradient[7]= 0.0;
//
//        gradient[8]= 0.0;
//        gradient[9]= 0.0;
//        gradient[10]= 1.0;
//        gradient[11]= -1.0;

        gradient[0]=  3 * a1 * (a1*x[0] + b1) * (a1*x[0] + b1);
        gradient[1]= -1.0;

        gradient[2]=  3 * a2 * (a2*x[0] + b2) * (a2*x[0] + b2);
        gradient[3]= -1.0;
    }
    result[0] = ((a1*x[0] + b1) * (a1*x[0] + b1) * (a1*x[0] + b1) - x[1]);
    result[1] = ((a2*x[0] + b2) * (a2*x[0] + b2) * (a2*x[0] + b2) - x[1]);
}

int main(int argc,char ** argv)
{
    ros::init(argc,argv,"motionprediction");
    ros::NodeHandle nh;
    ros::Rate rate(100);
    ros::Publisher predictionarray_pub=nh.advertise<offb_posctl::wallstate>("predictedwall_array",1, true);
    ros::Publisher predictionpath_pub=nh.advertise<nav_msgs::Path>("predictedwall_path",1);

    int timenodenumber=2;
    int constraintsnumber=3;
    double T=0.6;
    int timeconsumption=0;
    nlopt::opt opt(nlopt::LD_SLSQP, timenodenumber);
//    nlopt::opt opt(nlopt::LD_MMA, 2);
    chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now();
//    std::vector<double> lb(timenodenumber);

    std::vector<double> lb(2);
    lb[0] = -HUGE_VAL; lb[1] = 0;

    std::vector<double> ub(timenodenumber);
    std::vector<double> constraints_tol(constraintsnumber);

    std::vector<double> constraints_toltest(2);
    constraints_toltest[0]=1e-8;
    constraints_toltest[1]=1e-8;

//    lb[0] =0; lb[1] = 0; lb[2] =0; lb[3] = 0;
    ub[0] =T; ub[1] = T; ub[2] =T; ub[3] = T;
    constraints_tol[0]=1e-8; constraints_tol[1]=1e-8; constraints_tol[2]=1e-8;
    opt.set_lower_bounds(lb);
//    opt.set_upper_bounds(ub);
    opt.set_min_objective(myvfunc, NULL);

//    opt.add_inequality_constraint(myvconstraint, &data[0], 1e-8);
//    opt.add_inequality_constraint(myvconstraint, &data[1], 1e-8);
    opt.add_inequality_mconstraint(multiconstraint,data,constraints_toltest);
    opt.set_xtol_rel(1e-4);
    std::vector<double> x(2);
    x[0] = 1.234; x[1] = 5.678;
    double minf;

    try{
        nlopt::result result = opt.optimize(x, minf);
        std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
                  << std::setprecision(10) << minf << std::endl;
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
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