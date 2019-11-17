// conf_compare.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

/*changelog
# 0.4
findNearestsIndex in comf_compare.cpp is implemented

# 0.3
* positional input filenames can represent lists if they end with ini, comparison between the two files lines by lines; if positional input filename end with dconf, single comparison happens
* output file path is added as optional argument, program uses cout if parameter is not present
* find-to-compare switch is included but not implemented yet

# 0.2
difference in y values are allowed up to individual tolerance

# 0.1
First release
*/

#pragma region header with functions

#define VERSION_conf_compare 0.3

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include <boost/program_options.hpp> // to read in program call arguments
#include <boost/program_options/options_description.hpp> // to add descriptions of the program call arguments

namespace bpo = boost::program_options;
using disl = std::tuple<double, double, int>; // a dislocation is a (double, double, int) tuple for (posx,posy,type)

// evaluate ifname and put one or more filenames to ifnames
std::vector<std::string> processInputFile(std::string ifname)
{
    std::vector<std::string> ifnames;
    if (ifname.size() >= 4 && ifname.compare(ifname.size() - 4, 4, ".ini") == 0)
    {
        std::ifstream ifile(ifname);
        if (!ifile)
        {
            std::cerr << "The program couldn't open " << ifname << " for reading. Program terminates." << std::endl;
            exit(-1);
        }

        for (std::string tmp; ifile >> tmp;)
        {
            if (tmp.size() >= 6 && tmp.compare(tmp.size() - 6, 6, ".dconf") == 0)
                ifnames.push_back(tmp);
            else
            {
                std::cerr << ifname << " contains invalid dislocation configuration filename " << tmp << ". Program termiantes." << std::endl;
                exit(2);
            }
        }
        std::cout << ifname << " as input path contained " << ifnames.size() << " number of files to read." << std::endl;
        if (ifnames.size() == 0)
        {
            std::cerr << "Not enough input files. Program terminates." << std::endl;
            exit(1);
        }
    }
    else if (ifname.size() >= 6 && ifname.compare(ifname.size() - 6, 6, ".dconf") == 0)
    {
        std::cout << ifname << " will be used to read dislocation configuration" << std::endl;
        ifnames.push_back(ifname); // ifnames contains exactly 1 file that can be opened or
    }
    else
    {
        std::cerr << ifname << " is not a valid input filename. Program terminates." << std::endl;
        exit(1);
    }

    return ifnames;
}

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

// returns the index from the vector for which selectFrom[index] is closest to val
size_t findNearestsIndex(double val, const std::vector<double>& selectFrom)
{
    double larger = INFINITY;   // the smallest larger value in the list
    size_t largerID = 0;        // the position of larger in the list
    double smaller = 0;         // the largest smaller value in the list
    size_t smallerID = 0;       // the position of smaller in the list
    for (size_t  i = 0; i < selectFrom.size(); ++i)
    {
        if (selectFrom[i] > smaller && selectFrom[i] <= val)
        {
            smaller = selectFrom[i];
            smallerID = i;
        }
        if (selectFrom[i] < larger && selectFrom[i] >= val)
        {
            larger = selectFrom[i];
            largerID = i;
        }

    }
    if (val - smaller < larger - val)
        return smallerID;
    else
        return largerID;
}

#pragma endregion

int main(int argc, char** argv)
{
    std::vector<std::string> ifnames; // to store the two positional program call arguments of the input filenames

#pragma region reading in variables
    bpo::options_description requiredOptions("Required options"); // must be set from command line or file
    bpo::options_description optionalOptions("Optional options"); // must be set from command line or file

    requiredOptions.add_options()
        ("input-files", bpo::value<std::vector<std::string>>(&ifnames), "The input files to compare. No switchs are needed, they are the 1st and 2nd positional arguments. Files must end either with dconf for singular comparison or ini containing list of files.");

    bpo::positional_options_description positionalOptions;
    positionalOptions.add("input-files", 2);

    optionalOptions.add_options()
        ("output-filename,O", bpo::value<std::string>(), "The path where the result of the comparison(s) should be stored. If the value is not present, standard output will be used.")
        ("find-to-compare,f", "If set, input filename must point to file lists and files from the 1st list will be compared from a file from the 2nd list. Filenames will be treated as floating point values x and each file from the 1st list will compared with a file from the 2nd list that has the closest value to x.")
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

#pragma region processing input variables
    auto ifnames_a = processInputFile(ifnames[0]);
    auto ifnames_b = processInputFile(ifnames[1]);
    std::vector<double> ifnames_a_values;
    std::vector<double> ifnames_b_values;

    bool findToCompare = false;
    if (vm.count("find-to-compare"))
    {
        findToCompare = true;
        ifnames_a_values.resize(ifnames_a.size());
        ifnames_b_values.resize(ifnames_b.size());
        std::transform(ifnames_a.begin(), ifnames_a.end(), ifnames_a_values.begin(), [](std::string fname)
            {
                fname.erase(fname.size() - 6, 6);   // removes .dconf
                return std::stod(fname);            // interpret the filename as double
            });
        std::transform(ifnames_b.begin(), ifnames_b.end(), ifnames_b_values.begin(), [](std::string fname)
            {
                fname.erase(fname.size() - 6, 6);
                return std::stod(fname);
            });
    }

    if (ifnames_a.size() != ifnames_b.size() && !findToCompare)
    {
        std::cerr << "Error! The 1st and 2nd input files contain different number of dislocation configurations (" << ifnames_a.size() << "!=" << ifnames_b.size() << ") to process. Program terminates." << std::endl;
        exit(-1);
    }

    size_t fileListSize = ifnames_a.size();

    std::ofstream ofile;
    std::streambuf* coutbuf = std::cout.rdbuf(); //save old buf
    if (vm.count("output-filename"))
    {
        std::string ofname = vm["output-filename"].as<std::string>();
        ofile.open(ofname);
        if (!ofile)
        {
            std::cerr << "Error! Cannot open output file " << ofname << ". Program terminates.";
            exit(-1);
        }
        std::cout.rdbuf(ofile.rdbuf());
    }

    if (fileListSize > 1)
        std::cout << "# filename_a\tfilename_b\tlargest difference\tat ID\tat y\taverage difference\taverage difference square" << std::endl;


#pragma endregion

    for (unsigned int i = 0; i < fileListSize; ++i)
    {
        std::vector<disl> dislocsA, dislocsB; // containers of the N number of dislocations

        std::string ifname_a = ifnames_a[i];
        std::string ifname_b = ifnames_b[i];
        if (findToCompare)
            ifname_b = ifnames_b[findNearestsIndex(ifnames_a_values[i], ifnames_b_values)];

        readInDislocs(ifname_a, dislocsA);
        readInDislocs(ifname_b, dislocsB);

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

        double sum_fabs = 0;    // the sum of the absolute value of the normalized differences of the x position
        double sum_fabsSQ = 0;  // the sum of the square of the normalized differences of the x position
        double max_diff = 0;    // the largest absolute-normalized-difference of the x position
        size_t max_ID = 0;      // the ID of the dislocation with the highest difference

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
            sum_fabsSQ += fabs_diff * fabs_diff;

            if (fabs_diff > max_diff)
            {
                max_diff = fabs_diff;
                max_ID = i;
            }
        }
        double avg_fabs = sum_fabs / size;
        double avg_fabsSQ = sum_fabsSQ / size;

        if (fileListSize == 1)
        {
            if (max_diff > 1e-16)
            {
                std::cout << "\tLargest x-coordinate difference is \n"
                    << "\t\t d = " << max_diff << " for dislocation with \n"
                    << "\t\tID = " << max_ID << " and y coordinate\n"
                    << "\t\t y = " << std::get<1>(dislocsA[max_ID]) << "\n";
            }

            if (max_diff < indTol)
                std::cout << "All dislocations are at the same place.\n";

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
        }
        else
        {
            std::cout
                << ifnames_a[i] << "\t"
                << ifnames_b[i] << "\t"
                << max_diff << "\t"
                << max_ID << "\t"
                << std::get<1>(dislocsA[max_ID]) << "\t"
                << avg_fabs << "\t"
                << avg_fabsSQ << "\n";
        }
    }
    
    std::cout.rdbuf(coutbuf);
    std::cout << fileListSize << " pair of files have been compared." << std::endl;
    return 0;
}