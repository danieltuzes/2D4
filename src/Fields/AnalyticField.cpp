/*
 * SDDDST Simple Discrete Dislocation Dynamics Toolkit
 * Copyright (C) 2015-2019  Gábor Péterffy <peterffy95@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

#include "Fields/AnalyticField.h"

using namespace sdddstCore;

// calculates the stress of a dislocation wall
double f(double dx, double cos2piy, double dy)
{
    if (dx * dx + dy * dy > 1e-6)
    {
        double cosh2pix = cosh(dx * 2 * M_PI);
        double coshpixyd = cosh2pix - cos2piy;  // cosh2pix - cos2piy
        return  dx * (cosh2pix * cos2piy - 1) / (coshpixyd * coshpixyd) * 2 * M_PI * M_PI;
    }

    double dx2pi = dx * 2 * M_PI;
    double dy2pi = dy * 2 * M_PI;
    double cosh2pixminus1 = pow(dx2pi, 6) / 720 + pow(dx2pi, 4) / 24 + dx2pi * dx2pi * 0.5;
    double cos2piyminus1 = -pow(dy2pi, 6) / 720 + pow(dy2pi, 4) / 24 - dy2pi * dy2pi * 0.5;
    double coshminuscos = (pow(dy2pi, 6) + pow(dx2pi, 6)) / 720 + (pow(dx2pi, 4) - pow(dy2pi, 4)) / 24 + (dx2pi * dx2pi + dy2pi * dy2pi) * 0.5;
    return ((cosh2pixminus1 * cos2piyminus1 + cos2piyminus1 + cosh2pixminus1) * dx) / (coshminuscos * coshminuscos * 2 * M_PI * M_PI);
}

// calculates the stress of a images number of dislocation walls
template<int images>
double g(double dx, double cos2piy, double dy)
{
    return g<images - 1>(dx, cos2piy, dy) + f(dx - images, cos2piy, dy) + f(dx + images, cos2piy, dy);
}

template<>
double g<0>(double dx, double cos2piy, double dy)
{
    return f(dx, cos2piy, dy);
}

double f_dx(double dx, double cos2piy, double dy)
{
    double cosh2pix = cosh(M_PI * 2 * dx);
    double sinh2pix = sinh(M_PI * 2 * dx);

    if (dx * dx + dy * dy > 1e-6)
    {
        double coshpixyd = cosh2pix - cos2piy;  // cosh2pix - cos2piy
        return ((cosh2pix * cos2piy - 1) / (coshpixyd * coshpixyd) +
            dx * (sinh2pix * cos2piy * 2 * M_PI / (coshpixyd * coshpixyd) -
            (cosh2pix * cos2piy - 1) / (coshpixyd * coshpixyd * coshpixyd) * 4 * M_PI * sinh2pix)) * M_PI * M_PI * 2;
    }

    double dx2pi = dx * 2 * M_PI;
    double dy2pi = dy * 2 * M_PI;
    double cosh2pixminus1 = pow(dx2pi, 6) / 720 + pow(dx2pi, 4) / 24 + pow(dx2pi, 2) * 0.5;
    double cos2piyminus1 = -pow(dy2pi, 6) / 720 + pow(dy2pi, 4) / 24 - pow(dy2pi, 2) * 0.5;
    double coshminuscos = (pow(dy2pi, 6) + pow(dx2pi, 6)) / 720 + (pow(dx2pi, 4) - pow(dy2pi, 4)) / 24 + (pow(dx2pi, 2) + pow(dy2pi, 2)) * 0.5;
    double cosysinhx = -(pow(dx2pi, 7) * pow(dy2pi, 6)) / 3628800
        - (pow(dx2pi, 5) * pow(dy2pi, 6)) / 86400
        - (pow(dx2pi, 3) * pow(dy2pi, 6)) / 4320
        - (pow(dy2pi, 6) * dx2pi) / 720
        + (pow(dx2pi, 7) * pow(dy2pi, 4)) / 120960
        + (pow(dx2pi, 5) * pow(dy2pi, 4)) / 2880
        + (pow(dx2pi, 3) * pow(dy2pi, 4)) / 144
        + (pow(dy2pi, 4) * dx2pi) / 24
        - (pow(dx2pi, 7) * pow(dy2pi, 2)) / 1080
        - (pow(dx2pi, 5) * pow(dy2pi, 2)) / 240
        - (pow(dx2pi, 3) * pow(dy2pi, 2)) / 12
        - (pow(dy2pi, 2) * dx2pi) * 0.5
        + pow(dx2pi, 7) / 5040
        + pow(dx2pi, 5) / 120
        + pow(dx2pi, 3) / 6
        + dx2pi;
    return (
        (cosh2pixminus1 * cos2piyminus1 + cos2piyminus1 + cosh2pixminus1) / (coshminuscos * coshminuscos) +
        ((cosysinhx * dx) / (coshminuscos * coshminuscos) -
        (cosh2pixminus1 * cos2piyminus1 + cos2piyminus1 + cosh2pixminus1) / (coshminuscos * coshminuscos * coshminuscos) * dx * sinh(dx2pi) * 2) * M_PI * 2
        ) * 2 * M_PI * M_PI;


}

template<int images>
double g_dx(double dx, double cos2piy, double dy)
{
    return g_dx<images - 1>(dx, cos2piy, dy) + f_dx(dx - images, cos2piy, dy) + f_dx(dx + images, cos2piy, dy);
}

template<>
double g_dx<0>(double dx, double cos2piy, double dy)
{
    return f_dx(dx, cos2piy, dy);
}


AnalyticField::AnalyticField() :
    Field()
{
    // Nothing to do
}

AnalyticField::~AnalyticField()
{
    // Nothing to do
}

double AnalyticField::xy(double dx, double dy)
{
    double cos2piy = cos(M_PI * 2 * dy);
    if (dx < 0)
    {
        return -dx * f(dx + ANALYTIC_FIELD_N + 1, cos2piy, dy) +
            (1 + dx) * f(dx - ANALYTIC_FIELD_N, cos2piy, dy) +
            f(dx + ANALYTIC_FIELD_N, cos2piy, dy) +
            g<ANALYTIC_FIELD_N - 1>(dx, cos2piy, dy);
    }

    return dx * f(dx - ANALYTIC_FIELD_N - 1, cos2piy, dy) +
        (1 - dx) * f(dx + ANALYTIC_FIELD_N, cos2piy, dy) +
        f(dx - ANALYTIC_FIELD_N, cos2piy, dy) +
        g<ANALYTIC_FIELD_N - 1>(dx, cos2piy, dy);
}

double AnalyticField::xy_diff_x(double dx, double dy)
{
    double cos2piy = cos(M_PI * 2 * dy);
    if (dx < 0)
    {
        return -f(dx + ANALYTIC_FIELD_N + 1, cos2piy, dy) -
            dx * f_dx(dx + ANALYTIC_FIELD_N + 1, cos2piy, dy) +
            f_dx(dx - ANALYTIC_FIELD_N, cos2piy, dy) * (1 + dx) +
            f(dx - ANALYTIC_FIELD_N, cos2piy, dy) +
            f_dx(dx + ANALYTIC_FIELD_N, cos2piy, dy) +
            g_dx<ANALYTIC_FIELD_N - 1>(dx, cos2piy, dy);
    }

    return  f(dx - ANALYTIC_FIELD_N - 1, cos2piy, dy) +
        dx * f_dx(dx - ANALYTIC_FIELD_N - 1, cos2piy, dy) +
        f_dx(dx + ANALYTIC_FIELD_N, cos2piy, dy) * (1 - dx) -
        f(dx + ANALYTIC_FIELD_N, cos2piy, dy) +
        f_dx(dx - ANALYTIC_FIELD_N, cos2piy, dy) +
        g_dx<ANALYTIC_FIELD_N - 1>(dx, cos2piy, dy);
}
