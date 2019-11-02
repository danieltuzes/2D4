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

#include "Fields/Field.h"

using namespace sdddstCore;

Field::Field() {}

#if USE_IEEE_HYPERBOLIC // calculating cosh and sinh IEEE compatible way, except -ffast_math

// calculates the stress of a dislocation wall
double f(double dx, double cos2piy, double dy)
{
    if (dx * dx + dy * dy > 1e-6) // dislocations with larger distance
    {
        double cosh2pix = cosh(dx * 2 * M_PI);
        double coshpixyd = cosh2pix - cos2piy;  // cosh2pix - cos2piy
        return  dx * (cosh2pix * cos2piy - 1) / std::pow(coshpixyd, 2) * 2 * M_PI * M_PI;
    }

    // closer dislocations, to avoid (1-1)/0 type singularity, but with a Taylor expansion, one gets 0/0 type singularity
    double dx2pi = dx * 2 * M_PI;
    double dy2pi = dy * 2 * M_PI;
    double cosh2pixminus1 = std::pow(dx2pi, 6) / 720 + std::pow(dx2pi, 4) / 24 + std::pow(dx2pi, 2) * 0.5;
    double cos2piyminus1 = -std::pow(dy2pi, 6) / 720 + std::pow(dy2pi, 4) / 24 - std::pow(dy2pi, 2) * 0.5;
    double coshminuscos = (std::pow(dy2pi, 6) + std::pow(dx2pi, 6)) / 720 + (std::pow(dx2pi, 4) - std::pow(dy2pi, 4)) / 24 + (std::pow(dx2pi, 2) + std::pow(dy2pi, 2)) * 0.5;
    return ((cosh2pixminus1 * cos2piyminus1 + cos2piyminus1 + cosh2pixminus1) * dx) / std::pow(coshminuscos, 2) * 2 * M_PI * M_PI;
}

// calculates the stress of images number of dislocation walls
template<int images>
double g(double dx, double cos2piy, double dy)
{
    return g<images - 1>(dx, cos2piy, dy) + f(dx - images, cos2piy, dy) + f(dx + images, cos2piy, dy);
}

// end of recursive template g
template<>
double g<0>(double dx, double cos2piy, double dy)
{
    return f(dx, cos2piy, dy);
}

// calculates the derivative of the stress of a dislocation wall
double f_dx(double dx, double cos2piy, double dy)
{
    double dx2pi = dx * 2 * M_PI;

    if (dx * dx + dy * dy > 1e-6)
    {
        double cosh2pix = cosh(dx2pi);
        double sinh2pix = sinh(dx2pi);
        double coshpixyd = cosh2pix - cos2piy;  // cosh2pix - cos2piy
        return ((cosh2pix * cos2piy - 1) / std::pow(coshpixyd, 2) +
            dx * (sinh2pix * cos2piy * 2 * M_PI / std::pow(coshpixyd, 2) -
            (cosh2pix * cos2piy - 1) / std::pow(coshpixyd, 3) * 4 * M_PI * sinh2pix)) * M_PI * 2 * M_PI;
    }

    // closer dislocations
    double dy2pi = dy * M_PI * 2;
    double cosh2pixminus1 = std::pow(dx2pi, 6) / 720 + std::pow(dx2pi, 4) / 24 + std::pow(dx2pi, 2) * 0.5;
    double cos2piyminus1 = -std::pow(dy2pi, 6) / 720 + std::pow(dy2pi, 4) / 24 - std::pow(dy2pi, 2) * 0.5;
    double coshminuscos = (std::pow(dy2pi, 6) + std::pow(dx2pi, 6)) / 720 + (std::pow(dx2pi, 4) - std::pow(dy2pi, 4)) / 24 + (std::pow(dx2pi, 2) + std::pow(dy2pi, 2)) * 0.5;
    double cosysinhx = -(std::pow(dx2pi, 7) * std::pow(dy2pi, 6)) / 3628800
        - (std::pow(dx2pi, 5) * std::pow(dy2pi, 6)) / 86400
        - (std::pow(dx2pi, 3) * std::pow(dy2pi, 6)) / 4320
        - (std::pow(dy2pi, 6) * dx2pi) / 720
        + (std::pow(dx2pi, 7) * std::pow(dy2pi, 4)) / 120960
        + (std::pow(dx2pi, 5) * std::pow(dy2pi, 4)) / 2880
        + (std::pow(dx2pi, 3) * std::pow(dy2pi, 4)) / 144
        + (std::pow(dy2pi, 4) * dx2pi) / 24
        - (std::pow(dx2pi, 7) * std::pow(dy2pi, 2)) / 1080
        - (std::pow(dx2pi, 5) * std::pow(dy2pi, 2)) / 240
        - (std::pow(dx2pi, 3) * std::pow(dy2pi, 2)) / 12
        - (std::pow(dy2pi, 2) * dx2pi) * 0.5
        + std::pow(dx2pi, 7) / 5040
        + std::pow(dx2pi, 5) / 120
        + std::pow(dx2pi, 3) / 6
        + dx2pi;
    return (
        (cosh2pixminus1 * cos2piyminus1 + cos2piyminus1 + cosh2pixminus1) / std::pow(coshminuscos, 2) +
        ((cosysinhx * dx) / std::pow(coshminuscos, 2) -
        (cosh2pixminus1 * cos2piyminus1 + cos2piyminus1 + cosh2pixminus1) / std::pow(coshminuscos, 3) * dx * sinh(dx2pi) * 2) * 2 * M_PI
        ) * M_PI * 2 * M_PI;
}

// calculates the derivative of the stress of a images number of dislocation walls
template<int images>
double g_dx(double dx, double cos2piy, double dy)
{
    return g_dx<images - 1>(dx, cos2piy, dy) + f_dx(dx - images, cos2piy, dy) + f_dx(dx + images, cos2piy, dy);
}

// end of recursive template g
template<>
double g_dx<0>(double dx, double cos2piy, double dy)
{
    return f_dx(dx, cos2piy, dy);
}

// calculates the force field with stric IEEE precision
double Field::xy(double dx, double dy) const
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

// drops exp2pix and cos2piy and calls xy(dx, dy)
double Field::xy(double dx, double dy, double exp2pix, double cos2piy) const
{
    return xy(dx, dy);
}

// calculates the derived field with stric IEEE precision
double Field::xy_diff_x(double dx, double dy) const
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

// drops exp2pix and cos2piy and calls xy_diff_x(dx, dy)
double Field::xy_diff_x(double dx, double dy, double exp2pix, double cos2piy) const
{
    return xy_diff_x(dx, dy);
}
#else // calculating cosh and sinh faster but less precise way

// calculates the stress of a dislocation wall in the first cell (0 for 0th image)
double f_0(double dx, double cos2piy, double dy, double cosh2pix)
{
    double dx2pi = M_PI * 2 * dx;

    if (dx * dx + dy * dy > 1e-6) // dislocations with larger distance
    {
        double coshpixyd = cosh2pix - cos2piy;  // cosh2pix - cos2piy
        return  dx2pi * (cosh2pix * cos2piy - 1) / std::pow(coshpixyd, 2) * M_PI;
    }

    // closer dislocations, to avoid (1-1)/0 type singularity, but with a Taylor expansion, one gets 0/0 type singularity
    double dy2pi = M_PI * 2 * dy;
    double cosh2pixminus1 = std::pow(dx2pi, 6) / 720 + std::pow(dx2pi, 4) / 24 + std::pow(dx2pi, 2) * 0.5;
    double cos2piyminus1 = -std::pow(dy2pi, 6) / 720 + std::pow(dy2pi, 4) / 24 - std::pow(dy2pi, 2) * 0.5;
    double coshminuscos = (std::pow(dy2pi, 6) + std::pow(dx2pi, 6)) / 720 + (std::pow(dx2pi, 4) - std::pow(dy2pi, 4)) / 24 + (std::pow(dx2pi, 2) + std::pow(dy2pi, 2)) * 0.5;
    return ((cosh2pixminus1 * cos2piyminus1 + cos2piyminus1 + cosh2pixminus1) * dx2pi) / std::pow(coshminuscos, 2) * M_PI;
}

// calculates the stress of a dislocation wall, distance is at least 1 (l for larger)
double f_l(double dx, double cos2piy, double cosh2pix)
{
    double dx2pi = M_PI * 2 * dx;
    double coshpixyd = cosh2pix - cos2piy;  // cosh2pix - cos2piy
    return dx2pi * (cosh2pix * cos2piy - 1) / std::pow(coshpixyd, 2) * M_PI;
}

// calculates the stress of 3 dislocation wall pairs
double g3(double dx, double cos2piy, double dy, double exp2pix)
{
    return
        f_l(dx + 3, cos2piy, cosh__2pi_xp3) + f_l(dx - 3, cos2piy, cosh__2pi_xm3) +
        f_l(dx + 2, cos2piy, cosh__2pi_xp2) + f_l(dx - 2, cos2piy, cosh__2pi_xm2) +
        f_l(dx + 1, cos2piy, cosh__2pi_xp1) + f_l(dx - 1, cos2piy, cosh__2pi_xm1) +
        f_0(dx + 0, cos2piy, dy, cosh__2pi_x__);
}

// calculates the derivative of the stress of a dislocation wall
double f_dx_0(double dx, double cos2piy, double dy, double cosh2pix, double sinh2pix)
{
    double dx2pi = dx * 2 * M_PI;

    if (dx * dx + dy * dy > 1e-6)
    {
        double coshpixyd = cosh2pix - cos2piy;  // cosh2pix - cos2piy
        return
            (
            (cosh2pix * cos2piy - 1) / std::pow(coshpixyd, 2) +
                dx2pi * (sinh2pix * cos2piy / std::pow(coshpixyd, 2) - (cosh2pix * cos2piy - 1) / std::pow(coshpixyd, 3) * 2 * sinh2pix)
                )
            * M_PI * 2 * M_PI;
    }

    // closer dislocations
    double dy2pi = dy * M_PI * 2;
    double cosh2pixminus1 = std::pow(dx2pi, 6) / 720 + std::pow(dx2pi, 4) / 24 + std::pow(dx2pi, 2) * 0.5;
    double cos2piyminus1 = -std::pow(dy2pi, 6) / 720 + std::pow(dy2pi, 4) / 24 - std::pow(dy2pi, 2) * 0.5;
    double coshminuscos = (std::pow(dy2pi, 6) + std::pow(dx2pi, 6)) / 720 + (std::pow(dx2pi, 4) - std::pow(dy2pi, 4)) / 24 + (std::pow(dx2pi, 2) + std::pow(dy2pi, 2)) * 0.5;
    double cosysinhx = -(std::pow(dx2pi, 7) * std::pow(dy2pi, 6)) / 3628800
        - (std::pow(dx2pi, 5) * std::pow(dy2pi, 6)) / 86400
        - (std::pow(dx2pi, 3) * std::pow(dy2pi, 6)) / 4320
        - (std::pow(dy2pi, 6) * dx2pi) / 720
        + (std::pow(dx2pi, 7) * std::pow(dy2pi, 4)) / 120960
        + (std::pow(dx2pi, 5) * std::pow(dy2pi, 4)) / 2880
        + (std::pow(dx2pi, 3) * std::pow(dy2pi, 4)) / 144
        + (std::pow(dy2pi, 4) * dx2pi) / 24
        - (std::pow(dx2pi, 7) * std::pow(dy2pi, 2)) / 1080
        - (std::pow(dx2pi, 5) * std::pow(dy2pi, 2)) / 240
        - (std::pow(dx2pi, 3) * std::pow(dy2pi, 2)) / 12
        - (std::pow(dy2pi, 2) * dx2pi) * 0.5
        + std::pow(dx2pi, 7) / 5040
        + std::pow(dx2pi, 5) / 120
        + std::pow(dx2pi, 3) / 6
        + dx2pi;
    return (
        (cosh2pixminus1 * cos2piyminus1 + cos2piyminus1 + cosh2pixminus1) / std::pow(coshminuscos, 2) +
        ((cosysinhx * dx) / std::pow(coshminuscos, 2) -
        (cosh2pixminus1 * cos2piyminus1 + cos2piyminus1 + cosh2pixminus1) / std::pow(coshminuscos, 3) * dx * sinh(dx2pi) * 2) * 2 * M_PI
        ) * M_PI * 2 * M_PI;
}

// calculates the derivative of the stress of a dislocation wall, no condition for closeness
double f_dx_l(double dx, double cos2piy, double cosh2pix, double sinh2pix)
{
    double dx2pi = dx * 2 * M_PI;

    double chx_cy = cosh2pix - cos2piy;       // cosh2pix - cos2piy
    double chxcy_1 = cosh2pix * cos2piy - 1;  // cosh2pix * cos2piy - 1
    return
        (
            chxcy_1 / std::pow(chx_cy, 2)
            +
            dx2pi * (sinh2pix * cos2piy / std::pow(chx_cy, 2) - chxcy_1 / std::pow(chx_cy, 3) * 2 * sinh2pix)
            )
        * M_PI * 2 * M_PI;
}

// calculates the derivative of the stress of a dislocation wall, no condition for closeness, cosh2pix = sinh2pix in this range
double f_dx_h(double dx, double cos2piy, double coshsinh2pix)
{
    double dx2pi = dx * 2 * M_PI;

    double chx_cy = coshsinh2pix - cos2piy;       // cosh2pix - cos2piy
    double chxcy_1 = coshsinh2pix * cos2piy - 1;  // cosh2pix * cos2piy - 1
    return
        (
            chxcy_1 / std::pow(chx_cy, 2)
            +
            dx2pi * (coshsinh2pix * cos2piy / std::pow(chx_cy, 2) - chxcy_1 / std::pow(chx_cy, 3) * 2 * coshsinh2pix)
            )
        * M_PI * 2 * M_PI;
}

// calculates the derivative of the stress of 3 dislocation walls
double g3_dx(double dx, double cos2piy, double dy, double exp2pix)
{
    return
        f_dx_0(dx + 0, cos2piy, dy, cosh__2pi_x__, sinh__2pi_x__) +
        f_dx_l(dx + 1, cos2piy, cosh__2pi_xp1, sinh__2pi_xp1) + f_dx_l(dx - 1, cos2piy, cosh__2pi_xm1, sinh__2pi_xm1) +
        f_dx_l(dx + 2, cos2piy, cosh__2pi_xp2, sinh__2pi_xp2) + f_dx_l(dx - 2, cos2piy, cosh__2pi_xm2, sinh__2pi_xm2) +
        f_dx_h(dx + 3, cos2piy, cosh__2pi_xp3) + f_dx_h(dx - 3, cos2piy, cosh__2pi_xm3);
}

// returns the field at point dx, dy knowing their hyperbolic and trigonometric function values
double Field::xy(double dx, double dy, double exp2pix, double cos2piy) const
{
    if (dx < 0)
    {
        return
            f_l(dx + 5, cos2piy, cosh__2pi_xp5) * (0 - dx) +
            f_l(dx - 4, cos2piy, cosh__2pi_xm4) * (1 + dx) +
            f_l(dx + 4, cos2piy, cosh__2pi_xp4) +
            g3(dx, cos2piy, dy, exp2pix);
    }

    return
        f_l(dx - 5, cos2piy, cosh__2pi_xm5) * (0 + dx) +
        f_l(dx + 4, cos2piy, cosh__2pi_xp4) * (1 - dx) +
        f_l(dx - 4, cos2piy, cosh__2pi_xm4) +
        g3(dx, cos2piy, dy, exp2pix);
}

// returns the derivative of the field in x direction at point dx, dy knowing their hyperbolicand trigonometric function values
double Field::xy_diff_x(double dx, double dy, double exp2pix, double cos2piy) const
{
    if (dx < 0)
    {
        return
            f_dx_h(dx + 5, cos2piy, cosh__2pi_xp5) * (0 - dx) + f_l(dx + 5, cos2piy, cosh__2pi_xp5) * (-1) +
            f_dx_h(dx - 4, cos2piy, cosh__2pi_xm4) * (1 + dx) + f_l(dx - 4, cos2piy, cosh__2pi_xm4) * (+1) +
            f_dx_h(dx + 4, cos2piy, cosh__2pi_xp4) +
            g3_dx(dx, cos2piy, dy, exp2pix);
    }

    return
        f_dx_h(dx - 5, cos2piy, cosh__2pi_xm5) * (0 + dx) + f_l(dx - 5, cos2piy, cosh__2pi_xm5) * (+1) +
        f_dx_h(dx + 4, cos2piy, cosh__2pi_xp4) * (1 - dx) + f_l(dx + 4, cos2piy, cosh__2pi_xp4) * (-1) +
        f_dx_h(dx - 4, cos2piy, cosh__2pi_xm4) +
        g3_dx(dx, cos2piy, dy, exp2pix);
}

// returns the field at point dx, dy without knowing their hyperbolic and trigonometric function values
double Field::xy(double dx, double dy) const
{
    double exp2pix = exp(M_PI * 2 * dx);
    double cos2piy = cos(M_PI * 2 * dy);
    return xy(dx, dy, exp2pix, cos2piy);
}

// returns the derivative of the field in x direction at point dx, dy without knowing their hyperbolic and trigonometric function values
double Field::xy_diff_x(double dx, double dy) const
{
    double exp2pix = exp(M_PI * 2 * dx);
    double cos2piy = cos(M_PI * 2 * dy);
    return xy_diff_x(dx, dy, exp2pix, cos2piy);
}
#endif