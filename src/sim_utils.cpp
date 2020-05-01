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

// condmat Eq. 18, the weight expressing the strongness of the implicit part
double weight(double subSum, char T)
{
    if (subSum > 0)                                     // if you modify this function, you must also modify the weightInv too
    {
        if (T == 'c')                                   // original version used by Peterffy, different from the written in the arxiv paper
            return std::pow(subSum / (subSum + 1), 2);

        if (T == 'p')                                   // the version written in the paper, different from the original version used by Peterffy
            return subSum / (subSum + 1);

        if (T == 'm')                                   // mathematically hinted, the original function is (1 - subSum + 1 / (1 + subSum) - 2 * exp(-subSum)) / (1 - subSum - 1 / (1 + subSum)), but approximated
        {
            double c = 0.412073;                        // best fit on the exact mathematical with the approximation function via c
            return 1 / (1 + 1.5 / (subSum + c * subSum * subSum));
        }
    }

    return 0;
}

// returns the (old) subSum from the (old) weight
double weightInv(double weight, char T)
{
    if (weight > 0)                                     // if you modify this function, you must also modify the weight too
    {
        if (T == 'c')
            return sqrt(weight) / (1 - sqrt(weight));   // original version used by Peterffy, different from the written in the arxiv paper

        if (T == 'p')
            return weight / (1 - weight);               // the version written in the paper, different from the original version used by Peterffy

        if (T == 'm')                                   // the original function is not invertable, instead, this function is inverted: 1 / (1 + 1.5 / (subSum + c * subSum^2)) using c = 0.412073, a fit on [0:5]
        {
            double c = 0.412073;
            return (1 - weight - sqrt(1 - 2 * weight + 6 * c * weight + weight * weight - 6 * c * weight * weight)) / (2 * (-c + c * weight));
        }
    }

    return 0;
}