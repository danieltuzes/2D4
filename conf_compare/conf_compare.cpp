// conf_compare.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

/*changelog
# 1.2
deep comparison

# 1.1
bugfix: findNearest didn't return the correct closest val if it was the 0th and smaller, but 0

# 1.0
bugfix: ifname_b was falsely selected from ifnames_b and not from the sorted ifnames_values_b if findToCompare

# 0.9
bugfix: max_IDsy was not calculated

# 0.8
if findToCompare, files with larger time values than maxVal are not compared

# 0.7
* bugfix: program didn't use the filename for filename_a from sorted ifnames_a_values
* restructured code: comparison is moved to compareDislocConfs, renamed variables
* added feature: if findToCompare, then time value for file a and b are printed to the output

# 0.6
bugfix: program wrote out not the filename_b it was working with

# 0.5
if findToCompare, filenames are allowed to contain slashes and underscores before the time value: chars before the last slash and underscore are removed

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

#define VERSION_conf_compare 1.2

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <utility>

#include <boost/program_options.hpp> // to read in program call arguments
#include <boost/program_options/options_description.hpp> // to add descriptions of the program call arguments

namespace bpo = boost::program_options;
using disl = std::tuple<double, double, int>; // a dislocation is a (double, double, int) tuple for (posx,posy,type)

// evaluate ifname and return a vector<strin> with one or more filenames
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
bool readInDislocs(std::string ifname, std::vector<disl>& dislocs)
{
    std::ifstream ifile(ifname);
    if (!ifile)
    {
        std::cerr << "Error: cannot open file " << ifname << "." << std::endl;
        return false;
    }

    double sum_b = 0; // the sum of the Burger's vector, it must be 0

    // Reading in dislocations from ifile
    for (double x; ifile >> x;)
    {
        double y, b;
        if (!(ifile >> y && ifile >> b))
        {
            std::cerr << "Error in " << ifname << ". Cannot read in the y coordinate and Burger's vector for an x coordinate with value ~ " << x << "." << std::endl;
            return false;
        }

        if (std::fabs(b - rint(b)) > 1e-5)
        {
            std::cerr << "Error in " << ifname << ". Burger's vector supposed to be an integer, -1 or 1, but value " << b << " is found in file." << std::endl;
            return false;
        }

        dislocs.emplace_back(x, y, static_cast<int>(b));
        sum_b += b;
    }

    if (sum_b)
    {
        std::cerr << "Error in " << ifname << ". The sum of the Burger's vector supposed to be 0, but it is " << sum_b << "." << std::endl;
        return false;
    }

    return true;
}

// transform the difference into the range [-0.5:0.5)
void normalize(double& n)
{
    while (n < -0.5)
        n += 1;

    while (n >= 0.5)
        n -= 1;
}

double dist(const disl& a, const disl& b)
{
    double distx = std::get<0>(a) - std::get<0>(b);
    double disty = std::get<1>(a) - std::get<1>(b);
    normalize(distx);
    normalize(disty);
    return sqrt(distx * distx + disty * disty);
}

// compares dislocation configurations in ifname_a and ifname_b, checks if y value difference is smaller than indTol, and put the result into max_diff (largest difference), max_ID (ID of the dislocation with max_diff), max_IDsy (it's y value), avg_fabs (average of the absolute value of the difference), avg_fabsSQ (average of the difference square). If sort, dislocations will be sorted before comparison
bool compareDislocConfs(std::string ifname_a, std::string ifname_b, double indTol, double& max_diff, size_t& max_ID, double& max_IDsy, double& avg_fabs, double& avg_fabsSQ, bool sort, bool deep)
{
    std::vector<disl> dislocsA, dislocsB; // containers of the N number of dislocations

    if (!readInDislocs(ifname_a, dislocsA))
    {
        std::cerr << "Error during reading in " << ifname_a << "." << std::endl;
        return false;
    }
    if (!readInDislocs(ifname_b, dislocsB))
    {
        std::cerr << "Error during reading in " << ifname_b << "." << std::endl;
        return false;
    }
    if (dislocsA.size() != dislocsB.size())
    {
        std::cerr << "Error: the number of dislocations are not equal." << std::endl;
        return false;
    }

    size_t size = dislocsA.size();

    if (sort)
    {
        std::sort(dislocsA.begin(), dislocsA.end(), [](const disl& a, const disl& b) {return (std::get<1>(a) + std::get<2>(a)) > (std::get<1>(b) + std::get<2>(b)); });
        std::sort(dislocsB.begin(), dislocsB.end(), [](const disl& a, const disl& b) {return (std::get<1>(a) + std::get<2>(a)) > (std::get<1>(b) + std::get<2>(b)); });
    }

    double sum_fabs = 0;    // the sum of the absolute value of the normalized differences of the x position
    double sum_fabsSQ = 0;  // the sum of the square of the normalized differences of the x position

    for (size_t i = 0; i < size; ++i)
    {
        double y_diff = std::fabs(std::get<1>(dislocsA[i]) - std::get<1>(dislocsB[i]));
        normalize(y_diff);  // handling periodic boundary conditions
        if (y_diff > indTol)
        {
            std::cerr << "Error: y values are different (" << std::get<1>(dislocsA[i]) << " != " << std::get<1>(dislocsB[i]) << ") for dislocation with id = " << i << " and x coordinate\n"
                << std::get<0>(dislocsA[i]) << " in " << ifname_a << " and\n"
                << std::get<0>(dislocsB[i]) << " in " << ifname_b << "." << std::endl;
            return false;
        }

        double fabs_diff = std::fabs(std::get<0>(dislocsA[i]) - std::get<0>(dislocsB[i]));
        normalize(fabs_diff);   // handling periodic boundary conditions
        sum_fabs += fabs_diff;
        sum_fabsSQ += fabs_diff * fabs_diff;

        if (fabs_diff > max_diff)
        {
            max_diff = fabs_diff;
            max_ID = i;
        }
    }

    avg_fabs = sum_fabs / size;
    avg_fabsSQ = sum_fabsSQ / size;
    max_IDsy = std::get<1>(dislocsA[max_ID]);

    {
        std::string ofname = ifname_a + "_VS_" + ifname_b;
        std::ofstream of(ifname_a + "_VS_" + ifname_b);
        if (!of)
        {
            std::cerr << "Cannot create file for deeper analysis " << ofname << std::endl;
            return true;
        }
        of << "# distance\tdist to nearest disl in a\tits ID\tdist to nearest disl in b\tits ID\n";
        for (size_t i = 0; i < size; ++i)
        {
            of << dist(dislocsA[i], dislocsB[i]) << "\t";
            double smallestDistA = INFINITY;
            double smallestDistB = INFINITY;
            size_t IDA = size + 1;
            size_t IDB = size + 1;
            for (size_t j = 0; j < size; ++j)
            {
                if (j != i)
                {
                    double tmpA = dist(dislocsA[i], dislocsA[j]);
                    double tmpB = dist(dislocsB[i], dislocsB[j]);
                    if (tmpA < smallestDistA)
                    {
                        smallestDistA = tmpA;
                        IDA = j;
                    }
                    if (tmpB < smallestDistB)
                    {
                        smallestDistB = tmpB;
                        IDB = j;
                    }
                }
            }
            of << smallestDistA << "\t" << IDA << "\t" << smallestDistB << "\t" << IDB << "\n";
        }
    }


    return true;
}

// returns the index and the value from the vector for which selectFrom[index] is closest to val
std::pair<size_t, double> findNearest(double val, const std::vector<double>& selectFrom)
{
    double larger = INFINITY;   // the smallest larger value in the list
    size_t largerID = 0;        // the position of larger in the list
    double smaller = -INFINITY; // the largest smaller value in the list
    size_t smallerID = 0;       // the position of smaller in the list
    for (size_t i = 0; i < selectFrom.size(); ++i)
    {
        if (selectFrom[i] > smaller&& selectFrom[i] <= val)
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
        return std::pair<size_t, double>(smallerID, smaller);
    else
        return std::pair<size_t, double>(largerID, larger);
}

// a filename - value pair
class fname_value
{
public:
    fname_value() : m_fname(""), m_value(0) {};
    fname_value(std::string fname) : m_fname(fname), m_value(deduceValue()) {};

    // for std::sort
    bool operator<(const fname_value& rhs)
    {
        return m_value < rhs.m_value;
    }

    double value() const
    {
        return m_value;
    }

    std::string fname() const
    {
        return m_fname;
    }

    // determines the value from the fname
    double deduceValue() const
    {
        auto f_cpy = m_fname;
        f_cpy.erase(m_fname.size() - 6, 6);         // removes .dconf

        size_t fnameS = f_cpy.find_last_of("/");    // fname without the directory, i.e, the base name, starts here
        if (fnameS != std::string::npos)
            f_cpy.erase(0, fnameS + 1);             // removes chars until base name

        size_t tvalS = f_cpy.find_last_of("_");     // where the time value starts
        if (tvalS != std::string::npos)
            f_cpy.erase(0, tvalS + 1);              // removes chars until base name

        double val = std::stod(f_cpy);          // interpret the filename as double
        return val;
    }
private:
    std::string m_fname;
    double m_value;
};
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
        ("sort,s", "The input dislocations need to be sorted by Burger's vector and y direction. If switch is not used, disloations must be in the same order in both files.")
        ("individual-tolerance", bpo::value<double>()->default_value(1e-8), "The absolute value of the difference below which two coordinates considered to be the same.")
        ("similarity-tolerance", bpo::value<double>()->default_value(1e-6), "The average absolute value of the differences below which two realisation are similar.");
    ("deep,d", "Deeper analysis: all distance difference will be printed out, and the nearest dislocation is also shown.");

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
    auto ifnames_a = processInputFile(ifnames[0]);  // contains all (1 or more of) the filenames for files that has to be compared to
    auto ifnames_b = processInputFile(ifnames[1]);  // contains all (1 or more of) the filenames for files that can be compared with
    std::vector<fname_value> ifnames_values_a;      // filenames and their corresponding time values for sorting; sorted
    std::vector<fname_value> ifnames_values_b;      // filenames and their corresponding time values to pass to ib_values for selecting, sorted
    std::vector<double> if_b_values;                // contains time values for b files for selecting; sorted

    bool findToCompare = false;
    double maxVal;                                  // the largest time value to make comparison if findToCompare

    if (vm.count("find-to-compare"))
    {
        findToCompare = true;
        for (auto fname : ifnames_a)
            ifnames_values_a.emplace_back(fname);
        for (auto fname : ifnames_b)
            ifnames_values_b.emplace_back(fname);
        std::sort(ifnames_values_a.begin(), ifnames_values_a.end());
        std::sort(ifnames_values_b.begin(), ifnames_values_b.end());
        maxVal = std::min(ifnames_values_a.back().value(), ifnames_values_b.back().value());

        for (auto val : ifnames_values_b)
            if_b_values.push_back(val.value());
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

#pragma endregion

    std::cerr.precision(16);
    std::cout.precision(16);

    std::cout << "# Program call: ";
    for (int i = 0; i < argc; ++i)
        std::cout << argv[i] << " ";
    std::cout << std::endl;

    if (fileListSize > 1)
    {
        std::cout << "# ifname_a\tifname_b\tlargest difference\tat ID\tat y\taverage difference\taverage difference square";
        if (findToCompare)
            std::cout << "\ttime_a\ttime_b";
        std::cout << std::endl;
    }

    int successfulread = 0;
    for (unsigned int i = 0; i < fileListSize; ++i)
    {
        std::string ifname_a = ifnames_a[i];
        std::string ifname_b = ifnames_b[i];

        double ifname_a_value = 0;
        double nearestB_value = 0;
        if (findToCompare)
        {
            auto next_ifname_values_a = ifnames_values_a[i];    // the upcoming ifname_a with its value
            ifname_a = next_ifname_values_a.fname();            // The upcoming ifname_a in the ordered list
            ifname_a_value = next_ifname_values_a.value();      // The time value of the upcoming ifname_a in the ordered list
            auto nearest_b = findNearest(ifname_a_value, if_b_values);
            size_t nearestIndex = nearest_b.first;              // the index for which the value from ib_values is closest to nextIfname_a.value()
            nearestB_value = nearest_b.second;                  // the value closest to nextIfname_a.value()
            ifname_b = ifnames_values_b[nearestIndex].fname();  // the best ifname_b
            if (ifname_a_value > maxVal)                        // exit if there is no reason to continue comparison
                break;
        }

        double max_diff = 0;    // the largest absolute-normalized-difference of the x position
        size_t max_ID = 0;      // the ID of the dislocation with the highest difference
        double max_IDsy = 0;    // the y value of the dislocation with largest position difference
        double avg_fabs = 0;    // the average of the absolute value of the normalized differences of the x position
        double avg_fabsSQ = 0;  // the average of the square of the normalized differences of the x position

        if (!compareDislocConfs(ifname_a, ifname_b, indTol, max_diff, max_ID, max_IDsy, avg_fabs, avg_fabsSQ, vm.count("sort"), vm.count("deep")))
        {
            std::cerr << "Error: cannot compare files " << ifname_a << " and " << ifname_b << ". This pair is skipped." << std::endl;
            continue;
        }
        successfulread++;

        if (fileListSize == 1)
        {
            if (max_diff > 1e-16)
            {
                std::cout << "\tLargest x-coordinate difference is \n"
                    << "\t\t d = " << max_diff << " for dislocation with \n"
                    << "\t\tID = " << max_ID << " and y coordinate\n"
                    << "\t\t y = " << max_IDsy << "\n";
            }

            if (max_diff < indTol)
                std::cout << "All dislocations are at the same place.\n";

            if (avg_fabs > 1e-16)
            {
                std::cout << "\tAverage position difference:\n"
                    << "\t\t" << avg_fabs << "\n";
                if (avg_fabs < vm["similarity-tolerance"].as<double>())
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
                << ifname_a << "\t"
                << ifname_b << "\t"
                << max_diff << "\t"
                << max_ID << "\t"
                << max_IDsy << "\t"
                << avg_fabs << "\t"
                << avg_fabsSQ;
            if (findToCompare)
                std::cout
                << "\t" << ifname_a_value
                << "\t" << nearestB_value;
            std::cout << "\n";
        }
    }

    std::cout.rdbuf(coutbuf);
    std::cout << successfulread << " pair of files have been compared." << std::endl;
    return 0;
}