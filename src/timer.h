//
//  Timer.h
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 12/20/23.
//
#pragma once

#include <chrono>
#include <iostream>
#include <nlohmann/json.hpp>
#include <fstream>

///profiling the time of the entire subdivision process.
///
///an array of 8 {total time,time spent on single function, time spent on double functions, time spent on triple functions time spent on double functions' zero crossing test, time spent on three functions' zero crossing test, total subdivision time, total splitting time}
const int timer_amount = 10;

/// saves the timing profile to a file
/// @param[in] time_labels            The labels of the timings.
/// @param[in] timings           The values of the timings, corresponding to the lables.
/// @param[in] filename            The name of the output file.
///
/// @return         Whether this saving procedure is successful.
bool save_timings(const std::string& filename,
                  const std::array<std::string, timer_amount>& time_label,
                  const std::array<double, timer_amount>& timings);

/// An array storing the most current time profling
extern std::array<double, timer_amount> profileTimer;

/// The labels for timing stats.
const std::array<std::string, timer_amount> time_label = {"total time: ",
    "get active multiples: ",
    "single func: ",
    "two func: ",
    "three func: ",
    "sub two func: ",
    "sub three func: ",
    "subdivision: ",
    "evaluations: ",
    "splitting: "
};

/// the enum for the timing labels.
enum timeProfileName{
    total_time,
    getActiveMuti,
    singleFunc,
    twoFunc,
    threeFunc,
    sub_twoFunc,
    sub_threeFunc,
    subdivision,
    evaluation,
    splitting
};

/// add the timer recorded to the profiling timer.
/// @param[in] profile          The most current time profile.
/// @param[in] timer            The time recorded from this temporary timer.
///
/// @return         The updated time profile. 
std::array<double, timer_amount> combine_timer (const std::array<double, timer_amount> &profile, const std::array<double, timer_amount> &timer);

template<typename Fn>
class Timer
{
public:
    Timer(timeProfileName name,
          Fn&& func
          )
    : m_Name(name), m_Func(func)
    {
        starterTime = std::chrono::high_resolution_clock::now();
    }
    
    ~Timer(){}
    
    void Stop()
    {
        auto stopperTime = std::chrono::high_resolution_clock::now();
        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(starterTime).time_since_epoch().count();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
        auto duration = end - start;
        double ms = duration * 0.001;
        switch (m_Name){
            case total_time:
                m_timeProfile[0] += ms;
                break;
            case getActiveMuti:
                m_timeProfile[1] += ms;
                break;
            case singleFunc:
                m_timeProfile[2] += ms;
                break;
            case twoFunc:
                m_timeProfile[3] += ms;
                break;
            case threeFunc:
                m_timeProfile[4] += ms;
                break;
            case sub_twoFunc:
                m_timeProfile[5] += ms;
                break;
            case sub_threeFunc:
                m_timeProfile[6] += ms;
                break;
            case subdivision:
                m_timeProfile[7] += ms;
                break;
            case evaluation:
                m_timeProfile[8] += ms;
                break;
            case splitting:
                m_timeProfile[9] += ms;
                break;
            default:
                std::cout << "no matching time profile identifier" << std::endl;
        }
        m_Func(m_timeProfile);
    }
    
    
private:
    timeProfileName m_Name;
    Fn m_Func;
    std::array<double, timer_amount> m_timeProfile = {0,0,0,0,0,0,0,0,0,0};
    std::chrono::time_point<std::chrono::high_resolution_clock> starterTime;
};


