//
// Created by zm on 19-4-25.
//

#include "DOB.h"
#include <vector>
#include <queue>


DOB::DOB(int doblen,float dobrate) {
    DOB_data_in = {0.0f, 0.0f, 0.0f};
    DOB_list.push_back(DOB_data_in);
    dob_length=doblen;
    dob_rate=dobrate;
}

void DOB::add_data(float curtime, float curvel, float desacc) {
    DOB_data_in.cur_time = curtime;
    DOB_data_in.cur_vel = curvel;
    DOB_data_in.des_acc = desacc;

    if(DOB_list.size() == 1){
        delta_time = curtime;
    } else{
        delta_time = curtime - DOB_list.rbegin()->cur_time;
    }

    if(DOB_list.size() < dob_length){
        DOB_list.push_back(DOB_data_in);
    } else{
        std::vector<DOB_DATA>::iterator k_beg = DOB_list.begin();
        DOB_list.erase(k_beg);
        DOB_list.push_back(DOB_data_in);
    }
}

float DOB::dob_output() {

    if(DOB_list.size() < dob_length){

        return 0.0f;
    }

    else{
        std::vector<DOB_DATA>::iterator dob_k;
        float acc_sum = 0.0f;
        for(dob_k = DOB_list.begin(); dob_k != DOB_list.end(); ++ dob_k){
            acc_sum = acc_sum + dob_k->des_acc;
        }
        float data2return = ( acc_sum * delta_time - (DOB_list.back().cur_vel - DOB_list.front().cur_vel) ) / (DOB_list.back().cur_time - DOB_list.front().cur_time);

        return data2return;
    }

}

float DOB::acc_dob_output(float curracc, float desacc) {
    double tempsum=0.0;
    dob_list.push_back(std::make_pair(curracc, desacc));
    if(DOB_list.size() > dob_length){
        std::vector<std::pair<double, double> > ::iterator dob_iter = dob_list.begin();
        dob_list.erase(dob_iter);
    }
    std::vector<std::pair<double, double> > ::iterator dob_iter = dob_list.begin();
    while(dob_iter!=dob_list.end())
    {
        tempsum=tempsum+((*dob_iter).second-(*dob_iter).first);
        dob_iter++;
    }
    return dob_rate*tempsum/dob_list.size();

}