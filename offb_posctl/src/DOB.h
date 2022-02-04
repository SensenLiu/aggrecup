//
// Created by zm on 19-4-25.
//

#ifndef OFFB_POSCTL_DOB_H
#define OFFB_POSCTL_DOB_H

#include <vector>


class DOB_DATA {
public:
    float cur_time;
    float cur_vel;
    float cur_acc;
    float des_acc;

};


class DOB {
    float delta_time;
    int dob_length=10;
    float dob_rate=0.1;
    DOB_DATA DOB_data_in;
    std::vector <DOB_DATA> DOB_list; // [1st time, 2nd vel, 3rd acc]
    std::vector <std::pair<double, double> > dob_list; //积分平均表


public:
    DOB(int doblen,float dobrate);
    void add_data(float curtime, float curvel, float desacc);
    float dob_output();
    float acc_dob_output(float curracc, float desacc);
};

#endif //OFFB_POSCTL_DOB_H
