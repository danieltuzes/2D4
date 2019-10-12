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
#include <fstream>
#include <limits>
#include <iterator>
#include <cmath>

#include <boost/program_options.hpp>

#include "utility.h"
#include "count.h"
#include "datafile.h"

namespace bpo = boost::program_options;

static void fill_counter(const std::string& filename, IntersectionCounter& counter);

int main(int argc, char *argv[]){
  bpo::options_description options("Options");
  options.add_options()
    ("timetable", bpo::value<std::string>(),
        "File with 'time dconf-file-name' lines")
    ("numpoints", bpo::value<IntersectionCounter::size_type>()->default_value(100),
        "Number of considered shift values")
    ("counterfile", bpo::value<std::string>(),
        "Optional file to write counter values")
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
  auto numPoints = vm["numpoints"].as<IntersectionCounter::size_type>();

  auto counter = IntersectionCounter(numPoints);

  fill_counter(filename, counter);

  if (vm.count("counterfile"))
  {
    auto filename = vm["counterfile"].as<std::string>();
    std::ofstream file(filename);
    if (!file) {
      std::cerr << "Error: failed to create counter file: " << filename << std::endl;
      throw std::runtime_error("failed to open counter file");
    }
    for (auto&& p : counter) {
      file << p;
    }
    if (!file) {
      std::cerr << "Error when writing counter file" << std::endl;
    }
  }

  return 0;
}

static void fill_counter(const std::string& filename, IntersectionCounter& counter)
{
  std::vector<double> lastX;
  int fileNum = 0;
  std::vector<double>::size_type lineNum;
  std::vector<double>::size_type dislCount;

  TimeTable table(filename);

  for (auto&& tr : table) {
    DConf dconf(tr.filename);
    lineNum = 0;
    for (auto&& dr: dconf) {
      if (fileNum != 0) {
        double last = lastX[lineNum];
        double curr = dr.x;
        counter.addInterval(Interval(last, curr));
      }
      lastX[lineNum] = dr.x;
      ++lineNum;
    }
    if (!dconf) {
      // lineNumber is 0 based, but we detect the error after the increment,
      // which fortunately gives us the normal 1 based line number
      std::cerr << "Error while reading dislocation configuration:" << tr.filename << " line " << lineNum << std::endl;
      throw std::runtime_error("read error in dconf");
    }
    if (fileNum == 0) {
      dislCount = lineNum;
    }
    else { // check if dislCount is the same in every file
      if (lineNum != dislCount) {
        std::cerr << "Error: dislocation count in " << tr.filename << " is different then in previous files" << std::endl;
        throw std::runtime_error("unexpected dislocation count");
      }
    }
    ++fileNum;
  }

  if (!table) {
    std::cerr << "Error while reading timetable " << filename << std::endl;
    throw std::runtime_error("read error in timetable");
  }
}

