// 
// sim_utils.h : utility function declarations for sim_utils.cpp, used by simulation.cpp; version and changelog info is contained in simulation.h

#ifndef SDDDST_CORE_SIM_UTILS_H
#define SDDDST_CORE_SIM_UTILS_H

#include <vector>
#include <numeric>
#include <chrono>

#pragma region old utilities
 // transforms x into the range [-0.5:0.5)
void normalize(double& n);

// calculates the absolute value square of a vector as in mathematics
double absvalsq(const std::vector<double>& input);

// returns wall time in seconds
double get_wall_time();

#pragma endregion


#endif