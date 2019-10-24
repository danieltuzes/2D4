// conf_compare.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#pragma region header

#define VERSION_conf_compare 0.2
/*changelog
0.2 
difference in y values are allowed up to individual tolerance

0.1
First release
*/

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include <boost/program_options.hpp> // to read in program call arguments
#include <boost/program_options/options_description.hpp> // to add descriptions of the program call arguments

namespace bpo = boost::program_options;
using disl = std::tuple<double, double, int>; // a dislocation is a (double, double, int) tuple for (posx,posy,type)

// reads in dislocations from file named ifname and checks for errors
void readInDislocs(std::string ifname, std::vector<disl>& dislocs)
{
    std::ifstream ifile(ifname);
    if (!ifile)
    {
        std::cerr << "Error: cannot open file " << ifname << ". Program terminates.\n";
        exit(-1);
    }

    double sum_b = 0; // the sum of the Burger's vector, it must be 0

    // Reading in dislocations from ifile
    for (double x; ifile >> x;)
    {
        double y, b;
        if (!(ifile >> y && ifile >> b))
        {
            std::cerr << "Error in " << ifname << ". Cannot read in the y coordinate and Burger's vector for an x coordinate with value ~ " << x << ". Program terminates." << std::endl;
            exit(-1);
        }

        if (fabs(b - rint(b)) > 1e-5)
        {
            std::cerr << "Error in " << ifname << ". Burger's vector supposed to be an integer, -1 or 1, but value " << b << " is found. Program terminates." << std::endl;
            exit(-1);
        }

        dislocs.emplace_back(x, y, static_cast<int>(b));
        sum_b += b;
    }

    if (sum_b)
    {
        std::cerr << "Error in " << ifname << ". The sum of the Burger's vector supposed to be 0, but it is " << sum_b << ". Program terminates." << std::endl;
        exit(-1);
    }
}

// transform the difference into the range [-0.5:0.5)
void normalize(double& n)
{
    while (n < -0.5)
        n += 1;

    while (n >= 0.5)
        n -= 1;
}

#pragma endregion

int main(int argc, char** argv)
{
    std::vector<std::string> ifnames;
#pragma region reading in variables
    bpo::options_description requiredOptions("Required options"); // must be set from command line or file
    bpo::options_description optionalOptions("Optional options"); // must be set from command line or file

    requiredOptions.add_options()
        ("input-files", bpo::value<std::vector<std::string>>(&ifnames), "The input files to compare. No switch are needed, they are the 1st and 2nd positional arguments.");

    bpo::positional_options_description positionalOptions;
    positionalOptions.add("input-files", 2);

    optionalOptions.add_options()
        ("sort,s", bpo::value<std::string>()->default_value("y"), "The input dislocations will be sorted by Burger's vector and y direction, u means unsorted.")
        ("individual-tolerance", bpo::value<double>()->default_value(1e-8), "The absolute value of the difference below which two coordinates considered to be the same.")
        ("similarity-tolerance", bpo::value<double>()->default_value(1e-6), "The average absolute value of the differences below which two realisation are similar.");

    bpo::options_description options; // the superior container of the options

    options.add(requiredOptions).add(optionalOptions).add_options()
        ("help", "show this help")
        ("hide-copyright,c", "hides the copyright notice from the standard output");

    bpo::variables_map vm; // the storage for the variables

    try
    {
        bpo::store(bpo::command_line_parser(argc, argv).options(options).positional(positionalOptions).run(), vm);
        if (vm.count("help")) // if the user is curious 
        {
            std::cout << options << std::endl;
            exit(0);
        }

        bpo::notify(vm);
    }
    catch (bpo::error & e)
    {
        std::cerr << e.what() << std::endl;
        exit(-1);
    }

    if (!vm.count("hide-copyright"))
    {
        std::cout << "conf_compare (version " << VERSION_conf_compare << ") from 2D4 - a 2D discrete dislocation dynamics simulation program toolset.\n"
            << "Copyright (C) Dániel Tüzes <tuzes@metal.elte.hu>\n";
    }

    if (ifnames.size() != 2)
    {
        std::cerr << "Exactly 2 filenames are expected. Program termiantes.\n";
        exit(-1);
    }

    double indTol = vm["individual-tolerance"].as<double>(); // individual tolerance
#pragma endregion

    std::vector<disl> dislocsA, dislocsB; // container of the N number of dislocations

    readInDislocs(ifnames[0], dislocsA);
    readInDislocs(ifnames[1], dislocsB);

    std::cerr.precision(16);
    std::cout.precision(16);

    if (dislocsA.size() != dislocsB.size())
    {
        std::cerr << "The number of dislocations in the two files are not equal. Program terminates.\n";
        exit(-1);
    }
    size_t size = dislocsA.size();

    if (!vm["sort"].as<std::string>().compare("y"))
    {
        std::sort(dislocsA.begin(), dislocsA.end(), [](const disl& a, const disl& b) {return (std::get<1>(a) + std::get<2>(a)) > (std::get<1>(b) + std::get<2>(b)); });
        std::sort(dislocsB.begin(), dislocsB.end(), [](const disl& a, const disl& b) {return (std::get<1>(a) + std::get<2>(a)) > (std::get<1>(b) + std::get<2>(b)); });
    }

    double sum_fabs = 0; // the sum of the absolute value of the normalized differences the Burger's vector
    double max_diff = 0; // the largest absolute-normalized-difference of the Burger's vector
    size_t outlierID = 0; // the ID of the dislocation with the highest difference

    for (size_t i = 0; i < size; ++i)
    {
        if (std::abs(std::get<1>(dislocsA[i]) - std::get<1>(dislocsB[i])) > indTol)
        {
            std::cerr << "y values are different (" << std::get<1>(dislocsA[i]) << " != " << std::get<1>(dislocsB[i]) << ") for dislocation with id = " << i << " and x coordinate\n"
                << std::get<0>(dislocsA[i]) << " in " << ifnames[0] << " and\n"
                << std::get<0>(dislocsB[i]) << " in " << ifnames[1] << ". Program terminates.\n";
            exit(-1);
        }

        double fabs_diff = std::fabs(std::get<0>(dislocsA[i]) - std::get<0>(dislocsB[i]));
        sum_fabs += fabs_diff;

        if (fabs_diff > max_diff)
        {
            max_diff = fabs_diff;
            outlierID = i;
        }
    }

    if (max_diff > 1e-16)
    {
        std::cout << "\tLargest x-coordinate difference is \n"
            << "\t\t d = " << max_diff << " for dislocation with \n"
            << "\t\tID = " << outlierID << " and y coordinate\n"
            << "\t\t y = " << std::get<1>(dislocsA[outlierID]) << "\n";
    }
    if (max_diff < indTol)
        std::cout << "All dislocations are at the same place.\n";

    double avg_fabs = sum_fabs / size;
    
    if (avg_fabs > 1e-16)
    {
        std::cout << "\tAverage position difference:\n"
            << "\t\t" << avg_fabs << "\n";
        if (sum_fabs / size < vm["similarity-tolerance"].as<double>())
            std::cout << "The configurations are similar.\n";
        else
            std::cout << "The configurations are not similar.\n";
    }
    else
        std::cout << "The configurations are identical.\n";

    

    return 0;
}