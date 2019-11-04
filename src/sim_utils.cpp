// 
// sim_utils.h : utility function declarations for sim_utils.cpp, used by simulation.cpp

#include "sim_utils.h"

void normalize(double& n)
{
    if (n < -0.5) // it shouldn't be smaller than -0.9999999, "while" could be changed to "if", but it isn't faster
        n += 1;

    if (n >= 0.5) // it shouldn't be larger than 0.99999999, "while" could be changed to "if", but it isn't faster
        n -= 1;
}

// calculates the absolute value square of a vector as in mathematics
double absvalsq(const std::vector<double>& input)
{
    return std::accumulate(input.begin(), input.end(), 0., [](double a, double b) {return a + b * b; });
}

// returns wall time in seconds
double get_wall_time()
{
    // From https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_start_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(t_start);
    auto t_start_se = t_start_ms.time_since_epoch();

    double time_in_ms = static_cast<double>(t_start_se.count());

    return time_in_ms / 1000.;
}