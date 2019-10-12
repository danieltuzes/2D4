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
#include <cmath>

#include "count.h"

Interval::Interval(double s, double e) :
  startpos(s),
  endpos(e)
{}

std::ostream& operator<<(std::ostream& st, IntersectSample& is) {
  st << is.pos << "\t" << is.intersectCount << std::endl;
  return st;
}

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

IntersectionCounter::size_type
IntersectionCounter::nextIndexInDir(size_type index, Dir dir) {
  size_type nextIndex = index;
  if (dir == POS) {
    if (nextIndex == priv_size)
      nextIndex = 0;
    else
      ++nextIndex;
  }
  else { // dir == NEG
    if (nextIndex == 0)
      nextIndex = priv_size - 1;
    else
      --nextIndex;
  }
  return nextIndex;
}

IntersectionCounter::size_type
IntersectionCounter::closestInDir(double pos, Dir dir) {
  double posInDistUnits = (pos - startpos) / distance; // between 0 and approximate priv_size
  double roundedPos;
  if (dir == POS)
    roundedPos = std::ceil(posInDistUnits);
  else
    roundedPos = std::floor(posInDistUnits);

  if (roundedPos >= priv_size)
    return 0;
  else
    return roundedPos;
}

void IntersectionCounter::addInterval(const Interval& interval) {
  double s = interval.startpos;
  double e = interval.endpos;
  // direction: +/-1
  Dir dir = (s<=e) ? POS : NEG;
  double halfrange = (endpos - startpos)/2.0;
  if (std::abs(e-s) > halfrange) {
    // closer in the other direction
    dir = (dir == POS) ? NEG : POS;
  }
  IntersectionCounter::size_type start_index, end_index;
  start_index = closestInDir(s, dir);
  end_index = closestInDir(e, dir);
  for (auto i = start_index; i != end_index; i = nextIndexInDir(i, dir)) {
    ++points[i].intersectCount;
  }
}
