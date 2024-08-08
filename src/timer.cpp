//
//  timer.cpp
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 6/20/24.
//

#include "timer.h"

std::array<double, timer_amount> combine_timer (const std::array<double, timer_amount> &profile, const std::array<double, timer_amount> &timer){
    std::array<double, timer_amount> ret;
    for (int i = 0; i < timer_amount; i++){
        ret[i] = profile[i] + timer[i];
    }
    return ret;
}


