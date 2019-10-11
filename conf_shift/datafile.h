/* Datafile handling routines.
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

#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <iterator>

struct Record {
};

struct DislRecord : Record {
  double x;
  double y;
  double b;
};

std::istream& operator>>(std::istream& st, DislRecord& dr) {
  st >> dr.x >> dr.y >> dr.b;
  return st;
}

struct TimeRecord : Record {
  double time;
  std::string filename;
};

std::istream& operator>>(std::istream& st, TimeRecord& tr) {
  st >> tr.time >> tr.filename;
  return st;
}

std::istream& operator>>(std::istream& st, TimeRecord& tr);

template <class Rec>
class Datafile {
public:
  typedef std::istream_iterator<Rec> iterator;

  Datafile(const std::string& filename) :
    stream(filename)
  {
    if (!stream) {
      std::cerr << "Error: failed to open datafile: " << filename << std::endl;
      throw std::runtime_error("failed to open datafile");
    }
  }
  
  iterator begin() {
    return iterator(stream);
  }
  
  iterator end() {
    return iterator();
  }
  
  bool eof() {
    return stream.eof();
  }
  
  bool good() {
    return stream.good();
  }
  
  bool bad() {
    return stream.bad();
  }
  
  bool fail() {
    return stream.fail();
  }

  bool operator!() {
    return !stream;
  }
  
  explicit operator bool() {
    return bool(stream);
  }

private:
  std::ifstream stream;
};


typedef Datafile<DislRecord> DConf;

typedef Datafile<TimeRecord> TimeTable;
