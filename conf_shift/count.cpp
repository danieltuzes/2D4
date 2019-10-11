/* Counting logic.
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

#include <iostream>
#include <vector>
#include <stdexcept>

#include "count.h"

static int bound(int v, int min, int max) {
  if (v < min)
    return min;
  else if (v > max)
    return max;
  else
    return v;
}

Interval::Interval(double s, double e) :
  startpos(s),
  endpos(e)
{}

IntersectionCounter::IntersectionCounter(size_type n, double s, double e) :
  points(n),
  priv_size(n)
{
  if (s <= e) {
    startpos = s;
    endpos = e;
  }
  else {
    startpos = e;
    endpos = s;
  }

  if (n <= 0) {
    std::cerr << "Error: size not positive " << n << std::endl;
    throw std::invalid_argument("size not positive");
  }
  else {
    distance = (endpos - startpos)/priv_size;
  }

  for (size_type i = 0; i < priv_size; ++i) {
    points[i].pos = startpos + i*distance;
    points[i].intersectCount = 0;
  }
}

IntersectionCounter::size_type IntersectionCounter::size() {
  return priv_size;
}

IntersectionCounter::reference IntersectionCounter::operator[](size_type index) {
  return points[index];
}

IntersectionCounter::iterator IntersectionCounter::begin() {
  return points.begin();
}

IntersectionCounter::iterator IntersectionCounter::end() {
  return points.end();
}

void IntersectionCounter::addInterval(const Interval& interval) {

}
