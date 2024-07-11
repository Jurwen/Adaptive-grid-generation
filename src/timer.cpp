//
//  timer.cpp
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 6/20/24.
//

#include "timer.h"

std::array<double, timer_amount> profileTimer = {0,0,0,0,0,0,0,0,0,0};

std::array<double, timer_amount> combine_timer (const std::array<double, timer_amount> &profile, const std::array<double, timer_amount> &timer){
    std::array<double, timer_amount> ret;
    for (int i = 0; i < timer_amount; i++){
        ret[i] = profile[i] + timer[i];
    }
    return ret;
}

bool save_timings(const std::string& filename,
                  const std::array<std::string, timer_amount>& time_label,
                  const std::array<double, timer_amount>& timings)
{
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str(),std::ios::app);
    //fout.open(filename.c_str(),std::ios::app);
    json jOut;
    for (size_t i = 0; i < timings.size(); ++i) {
        jOut[time_label[i]] = timings[i];
    }
    //
    fout << jOut << std::endl;
    fout.close();
    return true;
    
    
}
