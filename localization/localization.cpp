// localization.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

/*changelog
# 1
both methods are implemented, the two different ways to measure localization

# 0.4
* utf8 encoding
* d was printed out mistakenly

# 0.3
* normalize is used to measure distance
* x_0 = P is normalized

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

#define VERSION_localization 1

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

    std::cout << "# filename\tCoM\tr\tP\teta" << std::endl;

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

        double alpha = atan2(Y, X);                 // the center of mass, origin, x axis angle

        double CoM = alpha / (2 * M_PI) + 1. / 4;   // the center of mass projected to the unit circle
        normalize(CoM);
        double r = sqrt(X * X + Y * Y);             // the distance from the origin, a kind of measure of localisation

        double d = 1;                               // the sum of the distances from the plane P to all dislocations
        double P = 0;                               // the selected plane for which d is minimal
        double freq = 1 / (1000. + x_vals.size());  // the frequency of the measuring points
        for (double Pnom = -0.5; Pnom < 0.5; Pnom += 1 / (1000. + x_vals.size()))
        {
            double dnom = std::accumulate(x_vals.begin(), x_vals.end(), 0., [Pnom](double sum, double x) {double dist = Pnom - x; normalize(dist); return sum + fabs(dist); }) / x_vals.size();;
            if (dnom < d)
            {
                d = dnom;
                P = Pnom;
            }
        }

        // fine search
        for (double Pnom = P - freq; Pnom < P + freq; Pnom += freq / 1000)
        {
            double dnom = std::accumulate(x_vals.begin(), x_vals.end(), 0., [Pnom](double sum, double x) {double dist = Pnom - x; normalize(dist); return sum + fabs(dist); }) / x_vals.size();;
            if (dnom < d)
            {
                d = dnom;
                P = Pnom;
            }
        }

        double eta = 1 - 4 * d;                     // another measure of localization

        std::cout << ifname << "\t" << CoM << "\t" << r << "\t" << P << "\t" << eta << "\n";
    }
    return 0;
}