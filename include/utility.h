/*
 * SDDDST Simple Discrete Dislocation Dynamics Toolkit
 * Copyright (C) 2015-2019  Gábor Péterffy <peterffy95@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

#ifndef SDDDST_CORE_UTILITY_H
#define SDDDST_CORE_UTILITY_H

#define _USE_MATH_DEFINES
#include <math.h>
#include <ctime>
#include <chrono>

inline void normalize(double& n)
{
    while (n < -0.5) // it shouldn't be smaller than -0.9999999, "while" could be changed to "if", but it isn't faster
        n += 1;

    while (n >= 0.5) // it shouldn't be larger than 0.99999999, "while" could be changed to "if", but it isn't faster
        n -= 1;
}

inline double X(double x)
{
    return sin(2 * M_PI * x) * 0.5 / M_PI;
}

inline double X2(double x)
{
    return (1 - cos(2 * M_PI * x)) * 0.5 / M_PI / M_PI;
}

inline double E(double x, double y, double K)
{
    return exp(-K * (X2(x) + X2(y)));
}

inline double X_dx(double x)
{
    return cos(2 * M_PI * x);
}

inline double X2_dx(double x)
{
    return sin(2 * M_PI * x) / M_PI;
}

inline double E_dx(double x, double y, double K)
{
    return -E(x, y, K) * K * X2_dx(x);
}


// From:
// https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows
double get_wall_time()
{
    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_start_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(t_start);
    auto t_start_se = t_start_ms.time_since_epoch();

    double time_in_ms = static_cast<double>(t_start_se.count());

    return time_in_ms / 1000.;
}

#endif
