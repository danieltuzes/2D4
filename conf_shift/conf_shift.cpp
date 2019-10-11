/* conf_shift calculates a nice cutting position for dislocation.
 *
 * Copyright (C) 2019 Vandrus Zoltán
 *
 * This file is part of conf_shift.
 *
 * conf_shift is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <string>
#include <iostream>
#include <limits>
#include <iterator>
#include <cmath>

#include <boost/program_options.hpp>

#include "utility.h"
#include "count.h"
#include "datafile.h"

namespace bpo = boost::program_options;

static void fill_bins(const std::string& filename, CountOnBins& counter);
static bool turns(double lastTurn, double lastX, double x);

// minimum difference in x position taken as change
const double X_MINDIFF = 1e-10;

static bool appr_less(double x, double y) {
  return (y-x) > X_MINDIFF;
}

enum Dir {
  POS,
  NEG,
  STILL
};

Dir get_dir(double x, double y) {
  if (appr_less(x, y))
    return POS;
  else if (appr_less(y, x))
    return NEG;
  else
    return STILL;
}

int main(int argc, char *argv[]){
  bpo::options_description options("Options");
  options.add_options()
    ("timetable", bpo::value<std::string>(),
        "File with 'time dconf-file-name' lines")
    ("bincount", bpo::value<CountOnBins::index_type>()->default_value(100),
        "Number of bins")
    ("help", "show this help")
  ;
  bpo::positional_options_description positional_options;
  positional_options.add("timetable", 1);

  bpo::variables_map vm; // the storage for the variables
  try {
    bpo::store(bpo::command_line_parser(argc, argv).
               options(options).positional(positional_options).run(), vm);
  }
  catch (bpo::error & e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  if (vm.count("help")) // if the user is curious
  {
    std::cout << options << std::endl;
    return 0;
  }

  if (!vm.count("timetable"))
  {
    std::cerr << "Error: timetable parameter required, use --help for command usage." << std::endl;
    return 1;
  }

  auto filename = vm["timetable"].as<std::string>();
  auto binCount = vm["bincount"].as<CountOnBins::index_type>();

  auto counter = CountOnBins(binCount);

  fill_bins(filename, counter);

  return 0;
}

static bool turns(double lastTurn, double lastX, double x) {
  Dir dir1 = get_dir(lastTurn, lastX);
  Dir dir2 = get_dir(lastX, x);
  if ((dir1 == POS && dir2 == NEG) ||
      (dir1 == NEG && dir2 == POS))
    return true;
  else
    return false;
}

static void fill_bins(const std::string& filename, CountOnBins& counter)
{
  std::vector<double> lastTurn;
  std::vector<double> lastX;
  int fileNum = 0;
  std::vector<double>::size_type lineNumber;
  
  TimeTable table(filename);
  
  for (auto&& tr : table) {
    DConf dconf(tr.filename);
    lineNumber = 0;
    for (auto&& dr: dconf) {
      if (fileNum == 0) {
        // reading first points as turnpoints
        lastTurn[lineNumber] = dr.x;
      }
      else if (fileNum == 1){
        lastX[lineNumber] = dr.x;
      }
      ++lineNumber; 
    }
    if (!dconf) {
      // lineNumber is 0 based, but we detect the error after the increment,
      // which fortunately gives us the normal 1 based line number
      std::cerr << "Error while reading dislocation configuration:" << tr.filename << " line " << lineNumber << std::endl;
    }
    ++fileNum;
  }
  
  if (!table) {
    std::cerr << "Error while reading timetable " << filename << std::endl;
    throw std::runtime_error("read error in timetable");
  }
}

