//
//  Timer.h
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 12/20/23.
//

#include <chrono>
#include <iostream>

//profiling the time of the entire subdivision process.
//an array of 4 {total time, time spent on double functions, time spent on triple functions, time spent on zero crossing test}
//currently, the zero crossing test is using linear programming based on Gurobi package.


enum timeProfileName{
    total,
    twoFunc,
    threeFunc,
    gurobi
};

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
    
    ~Timer()
    {
        Stop();
    }
    
    void Stop()
    {
        auto stopperTime = std::chrono::high_resolution_clock::now();
        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(starterTime).time_since_epoch().count();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
        auto duration = end - start;
        double ms = duration * 0.001;
        //std::cout << ms << "ms" << std::endl;
        switch (m_Name){
            case total:
                m_timeProfile[0] += ms;
                break;
            case twoFunc:
                m_timeProfile[1] += ms;
                break;
            case threeFunc:
                m_timeProfile[2] += ms;
                break;
            case gurobi:
                m_timeProfile[3] += ms;
                break;
            default:
                std::cout << "no matching time profile identifier" << std::endl;
        }
        m_Func(m_timeProfile);
    }
    
    
private:
    timeProfileName m_Name;
    Fn m_Func;
    std::valarray<double> m_timeProfile = {0,0,0,0};
    std::chrono::time_point<std::chrono::high_resolution_clock> starterTime;
};



