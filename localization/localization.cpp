// localization.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

/*changelog
# 0.3
normalize is used to measure distance

# 0.2
abs was mistakenly used isntead of fabs, again, really?

# 0.1
* first version working with the initial idea, seems working

*/


#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#define _USE_MATH_DEFINES
#include <math.h>

#include <boost/program_options.hpp> // to read in program call arguments
#include <boost/program_options/options_description.hpp> // to add descriptions of the program call arguments
namespace bpo = boost::program_options;

#define VERSION_localization 0.3

void normalize(double& n)
{
    while (n < -0.5) // bad predictions can lead to values like n=-10'000 or so and using hyperbolic functions on them is problematic
        n += 1;

    while (n >= 0.5) // in rare cases, if is not enough, and in most cases, this is faster than fprem1: https://stackoverflow.com/questions/58803438/best-way-calculating-remainder-on-floating-points
        n -= 1;
}

int main(int argc, char** argv)
{
    std::vector<std::string> ifnames; // to store the two positional program call arguments of the input filenames

#pragma region reading in variables
    bpo::options_description requiredOptions("Required options"); // must be set from command line or file

    requiredOptions.add_options()
        ("input-files", bpo::value<std::vector<std::string>>(&ifnames), "The input files for which the localization and its plane will be determined and printed out. No switchs are needed.");

    bpo::positional_options_description positionalOptions;
    positionalOptions.add("input-files", -1);

    bpo::options_description options; // the superior container of the options

    options.add(requiredOptions).add_options()
        ("help", "show this help");

    bpo::variables_map vm; // the storage for the variables

    try
    {
        bpo::store(bpo::command_line_parser(argc, argv).options(options).positional(positionalOptions).run(), vm);
        if (vm.count("help")) // if the user is curious 
        {
            std::cout << "localization (version " << VERSION_localization << ") from 2D4 - a 2D discrete dislocation dynamics simulation program toolset.\n"
                << "Copyright (C) Dániel Tüzes <tuzes@metal.elte.hu>\n";
            std::cout << options << std::endl;
            exit(0);
        }

        bpo::notify(vm);
    }
    catch (bpo::error& e)
    {
        std::cerr << e.what() << std::endl;
        exit(-1);
    }

    if (ifnames.size() < 1)
    {
        std::cerr << "At least 1 filename is expected. Program termiantes.\n";
        exit(-1);
    }

#pragma endregion

    std::cout << "# filename\tx_0\teta" << std::endl;

    for (auto ifname : ifnames)
    {
        std::ifstream i(ifname);
        if (!i)
        {
            std::cerr << "# Warning: cannot open " << ifname << ". Program skips this file." << std::endl;
            continue;
        }
        std::vector<double> x_vals;
        for (double tmp, x; i >> x && i >> tmp && i >> tmp; )
            x_vals.push_back(x);

        double X;                                   // the center of weight for the x coordinates distributed on the unit circle
        double Y;                                   // the center of weight for the y coordinates distributed on the unit circle
        X = std::accumulate(x_vals.begin(), x_vals.end(), 0., [](double sum, double x) {return sum + sin(2 * M_PI * x); }) / x_vals.size();
        Y = std::accumulate(x_vals.begin(), x_vals.end(), 0., [](double sum, double x) {return sum - cos(2 * M_PI * x); }) / x_vals.size();

        //double a = sqrt(X * X + Y * Y);           // the distance of the center of mass from the unit circle's circumference, also a type of measure for localization
        double alpha = atan2(Y, X);                 // the center of mass, origin, x axis angle

        double P = (alpha + M_PI / 2) / (2 * M_PI); // the center of mass projected to the unit circle
        normalize(P);
        double d;                                   // the sum of the distances from the plane P to all dislocations
        d = std::accumulate(x_vals.begin(), x_vals.end(), 0., [P](double sum, double x) {double dist = P - x; normalize(dist); return sum + fabs(dist); }) / x_vals.size();

        double eta = 1 - 4 * d;

        std::cout << ifname << "\t" << P << "\t" << eta << "\t" << d << "\n";
    }
    return 0;
}