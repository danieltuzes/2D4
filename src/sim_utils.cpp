// 
// sim_utils.h : utility function declarations for sim_utils.cpp, used by simulation.cpp

#include "sim_utils.h"

void normalize(double& n)
{
    while (n < -0.5) // bad predictions can lead to values like n=-10'000 or so and using hyperbolic functions on them is problematic
        n += 1;

    while (n >= 0.5) // in rare cases, if is not enough, and in most cases, this is faster than fprem1: https://stackoverflow.com/questions/58803438/best-way-calculating-remainder-on-floating-points
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