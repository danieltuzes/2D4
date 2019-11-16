// merge_data_func.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

/*
# 0.1
First version tracket file
*/

#pragma region header namespaces, includes and defines
#define VERSION_merge_data_func 0.1

#include <iostream>
#include <fstream>

#include <boost/program_options.hpp> // to read in program call arguments
#include <boost/program_options/options_description.hpp> // to add descriptions of the program call arguments

namespace bpo = boost::program_options;

#pragma endregion

#pragma region functions and classes

// datapoints consisting a Tx x and an Ty y value; they can be compared based on their x value
template <typename Tx, typename Ty>
class DP
{
public:
    DP() : m_x(0), m_y() {};
    DP(Tx x, Ty y) : m_x(x), m_y(y) {};

    Tx x() const
    {
        return m_x;
    }

    Ty y() const
    {
        return m_y;
    }

    Tx& x()
    {
        return m_x;
    }

    Ty& y()
    {
        return m_y;
    }

    Tx& x(Tx x)
    {
        m_x = x;
        return m_x;
    }

    Ty& y(Ty y)
    {
        m_y = y;
        return m_y;
    }

    friend bool operator< (const DP<Tx, Ty>& lhs, const DP<Tx, Ty>& rhs)
    {
        return lhs.m_x < rhs.m_x;
    }
    friend bool operator> (const DP<Tx, Ty>& lhs, const DP<Tx, Ty>& rhs)
    {
        return lhs.m_x > rhs.m_x;
    }
    friend bool operator== (const DP<Tx, Ty>& lhs, const DP<Tx, Ty>& rhs)
    {
        return lhs.m_x == rhs.m_x;
    }
    friend bool operator<= (const DP<Tx, Ty>& lhs, const DP<Tx, Ty>& rhs)
    {
        return lhs.m_x <= rhs.m_x;
    }
    friend bool operator>= (const DP<Tx, Ty>& lhs, const DP<Tx, Ty>& rhs)
    {
        return lhs.m_x >= rhs.m_x;
    }

private:
    Tx m_x;
    Ty m_y;
};

// reads in the file fname and keep the x_col-th column as x and the y_cols as y values
std::vector<DP<double, std::vector<double>>> readFromFile(std::string fname, int x_col, std::vector<int> y_cols)
{
    std::ifstream ifile(fname);
    if (!ifile)
    {
        std::cerr << "Cannot open " << fname << " for reading. Program terminates.\n";
        exit(-1);
    }
    std::vector<DP<double, std::vector<double>>> data;

    for (std::string line; std::getline(ifile, line);)
    {
        if (line == "")
            continue;

        std::stringstream ss(line);

        double x = 0;                                       // the x value in that line
        std::vector<double> y_col_val(y_cols.size(),0);     // the different type of values for y in that line

        size_t lastYCol = 0;
        for (int i = 0; ; ++i)                              // iterating through the columns
        {
            std::string value;                              // a value of the actual column
            ss >> value;


            if (i == x_col)
                x = std::stod(value);

            auto y_col_match = std::find(y_cols.begin(), y_cols.end(), i); // the position of i in the list
            if (y_col_match != y_cols.end())                // if i is in the list
            {
                y_col_val[lastYCol] = std::stod(value);     // then store the column
                lastYCol++;
            }

            if (!ss)                                        // if cannot read more from that line
                break;
        }
        data.emplace_back(x, y_col_val);
    }
    return data;
}

#pragma endregion


int main(int argc, char** argv)
{
#pragma region reading in variables
    std::vector<std::string> ifnames;   // will store all the fnames, but at the beginning, may contain only the filename of the ini file
    std::vector<int> y_cols;            // the y values for which the averaging should be made

    bpo::options_description requiredOptions("Required options"); // must be set from command line or file
    bpo::options_description optionalOptions("Optional options"); // must be set from command line or file

    requiredOptions.add_options()
        ("input-files", bpo::value<std::vector<std::string>>(&ifnames)->composing(), "The path to an .ini file that contains the path to the files to merge, or the list of files themselves. No switch is needed, this is the first positional argument.");

    bpo::positional_options_description positionalOptions;
    positionalOptions.add("input-files", -1);

    optionalOptions.add_options()
        ("x-column,x", bpo::value<int>()->default_value(0), "The column index (starts from 0) for the x values.")
        ("y-columns,y", bpo::value<std::vector<int>>(&y_cols)->multitoken()->default_value(std::vector<int>{1}, "1"), "The column index(es) (start from 0) for the y values.")
        ("output-filename,O", bpo::value<std::string>(), "The path where the result of the comparison(s) should be stored. If the value is not present, standard output will be used.");

    bpo::options_description options;   // the superior container of the options

    options.add(requiredOptions).add(optionalOptions).add_options()
        ("help", "show this help")
        ("hide-copyright,c", "hides the copyright notice from the standard output");

    bpo::variables_map vm;              // the storage for the variables

    try
    {
        bpo::store(bpo::command_line_parser(argc, argv).options(options).positional(positionalOptions).run(), vm);
        if (vm.count("help"))           // if the user is curious 
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
        std::cout
            << "merge_data_func (version " << VERSION_merge_data_func << ") from 2D4 - a 2D discrete dislocation dynamics simulation program toolset.\n"
            << "Copyright (C) Dániel Tüzes <tuzes@metal.elte.hu>\n";
    }

    if (ifnames.empty())
    {
        std::cout << "Type in the name of the files to compare: \n";
        std::string ifnamelist;
        std::getline(std::cin, ifnamelist);
        std::stringstream ss(ifnamelist);
        for (std::string fname; ss >> fname; ifnames.push_back(fname));
        if (ifnames.empty())
        {
            std::cerr << "At least an ini file containing a list files or two filenames to compare should be provided. Program temrinates.\n";
            exit(-1);
        }
    }

#pragma endregion

#pragma region processing input variables
    
    // handling output file: if present, cout will print to that file

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

    // reading in ini file's input filenames to ifnames; ifname[0] will be the first input datafile's name
    if (ifnames.size() == 1)
    {
        if (ifnames[0].substr(ifnames[0].size() - 4, 4) != ".ini")
        {
            std::cerr << "If only 1 input file is given, it must be an ini file containing the list of files to compare. Program terminates.\n";
            exit(-1);
        }
        std::ifstream ifile(ifnames[0]);
        if (!ifile)
        {
            std::cerr << "Cannot open file " << ifnames[0] << " for reading. Program terminates.\n";
            exit(-1);
        }
        ifnames.clear();
        for (std::string line; ifile >> line; ifnames.push_back(line));
    }

    int x_col = vm["x-column"].as<int>();

#pragma endregion

    size_t files_to_merge = ifnames.size();                 // the number of files to merge
    std::vector<DP<double, std::vector<double>>> merged;    // the union of the data points
    std::vector<std::vector<DP<double, std::vector<double>>>> ifile_datas;

    // print out program call
    {
        std::cout << "# Program is called as: ";
        for (int i = 0; i < argc; ++i)
            std::cout << argv[i] << " ";
        std::cout << std::endl;
    }
    
    for (const auto fname : ifnames)
    {
        ifile_datas.push_back(readFromFile(fname, x_col, y_cols));
        for (const auto& xy : ifile_datas.back())
            merged.emplace_back(xy.x(), std::vector<double>(y_cols.size(), 0));
    }

    std::sort(merged.begin(), merged.end());            // calls friend bool operator< (const DP<Tx, Ty>& lhs, const DP<Tx, Ty>& rhs)
    auto it = std::unique(merged.begin(), merged.end());// remove duplicated elements
    merged.resize(it - merged.begin());                 // shrink to the unique elements
    std::vector<int> merged_c(merged.size(),0);         // counts how many values are added to the corresponding merged value

    for (const auto& ifile_data : ifile_datas)          // for each data file
    {
        int m_i = 0;                                    // index moving through the merged values
        for (size_t i = 0; i < ifile_data.size(); ++i)  // for each datapoint in the file
        {
            while (m_i < merged.size() && merged[m_i] < ifile_data[i])                  // skip the merged values till a datapoint is available
                m_i++;

            while (m_i < merged.size() && i < ifile_data.size()-1 && merged[m_i] <= ifile_data[i + 1])
            {
                double deltax = ifile_data[i + 1].x() - ifile_data[i].x();              // the distance in the x values for the two points on the left and right side of the merged value
                double weight_r = 1 - (ifile_data[i+1].x() - merged[m_i].x()) / deltax; // the weight of the right value from the linear interpolation
                double weight_l = 1 - (merged[m_i].x() - ifile_data[i].x()) / deltax;   // the weight of the left value from the linear interpolation

                for (int j = 0; j < ifile_data[i].y().size(); ++j)                      // interpolate all the requested y values
                    merged[m_i].y()[j] += weight_r * ifile_data[i+1].y()[j] + weight_l * ifile_data[i].y()[j];

                merged_c[m_i]++;                        // at averaging it will be needed to know with that it should be divided
                m_i++;                                  // next merged item, please
            }
        }
    }

    // divide the y values with the number of points added 
    for (size_t i = 0; i < merged.size(); ++i)
    {
        std::cout << merged[i].x();
        for (size_t j = 0; j < merged[i].y().size(); ++j)
            std::cout << "\t" << merged[i].y()[j] / merged_c[i];
        std::cout << std::endl;
    }

    // prints out to the real stdout
    std::cout.rdbuf(coutbuf);
    std::cout << files_to_merge << " files have been merged.\n";
    return 0;
}