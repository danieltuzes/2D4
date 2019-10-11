/* Counting logic header file.
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

#include <vector>

class Interval {
public:
  double startpos;
  double endpos;
  Interval(double s, double e);
};

struct IntersectSample {
  double pos;
  int intersectCount;
};

class IntersectionCounter {
private:
  typedef std::vector<IntersectSample> container_type;
public:
  typedef container_type::size_type size_type;
  typedef container_type::iterator iterator;
  typedef container_type::reference reference;
  // Count intersections on n different values s, ... s+(n-1)*distance
  IntersectionCounter(size_type n, double s=-0.5, double e=0.5);
  size_type size();
  reference operator[](size_type index);
  iterator begin();
  iterator end();
  void addInterval(const Interval& interval);
private:
  container_type points;
  size_type priv_size;
  // startpos <= endpos always
  double startpos;
  double endpos;
  double distance;
};
