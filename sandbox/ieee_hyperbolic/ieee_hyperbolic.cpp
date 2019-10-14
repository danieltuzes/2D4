// ieee_hyperbolic.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#define ANALYTIC_FIELD_N 4

#pragma region IEEE hyperbolic
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


double IEEE_xy(double dx, double dy)
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

double IEEE_xy_diff_x(double dx, double dy)
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
#pragma endregion

#pragma region hyperbolic with identities

#pragma region defines
#define cosh2pi 267.746761483748222245931879901    // cosh(2 * pi)
#define cosh4pi 143375.656570070321051450310133    // cosh(4 * pi)
#define cosh6pi 7.67764676977233502175191858009e7  // cosh(6 * pi)
#define cosh8pi 4.11131577927974976374895337541e10 // cosh(8 * pi)
#define coshTpi 2.20157529303160145057002722946e13 // cosh(10* pi)

#define sinh2pi 267.744894041016514257117449688    // sinh(2 * pi)
#define sinh4pi 143375.656566582978695241314640    // sinh(4 * pi)
#define sinh6pi 7.67764676977233437051070497210e7  // sinh(6 * pi), basically cosh6pi
#define sinh8pi 4.11131577927974976374773721974e10 // sinh(8 * pi), basically cosh8pi
#define sinhTpi 2.20157529303160145057002722719e13 // sinh(10* pi), basically coshTpi

#define cosh__2pi_xp1 (cosh2pix * cosh2pi + sinh2pix * sinh2pi)  // == cosh( 2*pi * (x+1) )
#define cosh__2pi_xm1 (cosh2pix * cosh2pi - sinh2pix * sinh2pi)  // == cosh( 2*pi * (x-1) )
#define cosh__2pi_xp2 (cosh2pix * cosh4pi + sinh2pix * sinh4pi)  // == cosh( 2*pi * (x+2) )
#define cosh__2pi_xm2 (cosh2pix * cosh4pi - sinh2pix * sinh4pi)  // == cosh( 2*pi * (x-2) )
#define cosh__2pi_xp3 (cosh2pix * cosh6pi + sinh2pix * sinh6pi)  // == cosh( 2*pi * (x+3) )
#define cosh__2pi_xm3 (cosh2pix * cosh6pi - sinh2pix * sinh6pi)  // == cosh( 2*pi * (x-3) )
#define cosh__2pi_xp4 (cosh2pix * cosh8pi + sinh2pix * sinh8pi)  // == cosh( 2*pi * (x+4) )
#define cosh__2pi_xm4 (cosh2pix * cosh8pi - sinh2pix * sinh8pi)  // == cosh( 2*pi * (x-4) )
#define cosh__2pi_xp5 (cosh2pix * coshTpi + sinh2pix * sinhTpi)  // == cosh( 2*pi * (x+5) )
#define cosh__2pi_xm5 (cosh2pix * coshTpi - sinh2pix * sinhTpi)  // == cosh( 2*pi * (x-5) )

#define sinh__2pi_xp1 (sinh2pix * cosh2pi + cosh2pix * sinh2pi)  // == sinh( 2*pi * (x+1) )
#define sinh__2pi_xm1 (sinh2pix * cosh2pi - cosh2pix * sinh2pi)  // == sinh( 2*pi * (x-1) )
#define sinh__2pi_xp2 (sinh2pix * cosh4pi + cosh2pix * sinh4pi)  // == sinh( 2*pi * (x+2) )
#define sinh__2pi_xm2 (sinh2pix * cosh4pi - cosh2pix * sinh4pi)  // == sinh( 2*pi * (x-2) )
#define sinh__2pi_xp3 (sinh2pix * cosh6pi + cosh2pix * sinh6pi)  // == sinh( 2*pi * (x+3) )
#define sinh__2pi_xm3 (sinh2pix * cosh6pi - cosh2pix * sinh6pi)  // == sinh( 2*pi * (x-3) )
#define sinh__2pi_xp4 (sinh2pix * cosh8pi + cosh2pix * sinh8pi)  // == sinh( 2*pi * (x+4) )
#define sinh__2pi_xm4 (sinh2pix * cosh8pi - cosh2pix * sinh8pi)  // == sinh( 2*pi * (x-4) )
#define sinh__2pi_xp5 (sinh2pix * coshTpi + cosh2pix * sinhTpi)  // == sinh( 2*pi * (x+5) )
#define sinh__2pi_xm5 (sinh2pix * coshTpi - cosh2pix * sinhTpi)  // == sinh( 2*pi * (x-5) )
#pragma endregion

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

// calculates the stress of a dislocation wall, distance is at least 1
double f_l(double dx, double cos2piy, double dy, double cosh2pix)
{
    double dx2pi = M_PI * 2 * dx;
    double coshpixyd = cosh2pix - cos2piy;  // cosh2pix - cos2piy
    return dx2pi * (cosh2pix * cos2piy - 1) / std::pow(coshpixyd, 2) * M_PI;
}

// calculates the stress of 3 dislocation wall pairs
double g3(double dx, double cos2piy, double dy, double cosh2pix, double sinh2pix)
{
    return
        f_0(dx + 0, cos2piy, dy, cosh2pix) +
        f_l(dx + 1, cos2piy, dy, cosh__2pi_xp1) + f_l(dx - 1, cos2piy, dy, cosh__2pi_xm1) +
        f_l(dx + 2, cos2piy, dy, cosh__2pi_xp2) + f_l(dx - 2, cos2piy, dy, cosh__2pi_xm2) +
        f_l(dx + 3, cos2piy, dy, cosh__2pi_xp3) + f_l(dx - 3, cos2piy, dy, cosh__2pi_xm3);
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
double f_dx_l(double dx, double cos2piy, double dy, double cosh2pix, double sinh2pix)
{
    double dx2pi = dx * 2 * M_PI;

    double coshpixyd = cosh2pix - cos2piy;  // cosh2pix - cos2piy
    return
        (
        (cosh2pix * cos2piy - 1) / std::pow(coshpixyd, 2) +
            dx2pi * (sinh2pix * cos2piy / std::pow(coshpixyd, 2) - (cosh2pix * cos2piy - 1) / std::pow(coshpixyd, 3) * 2 * sinh2pix)
            )
        * M_PI * 2 * M_PI;
}

// calculates the derivative of the stress of 3 dislocation walls
double g3_dx(double dx, double cos2piy, double dy, double cosh2pix, double sinh2pix)
{
    return
        f_dx_0(dx + 0, cos2piy, dy, cosh2pix, sinh2pix) +
        f_dx_l(dx + 1, cos2piy, dy, cosh__2pi_xp1, sinh__2pi_xp1) + f_dx_l(dx - 1, cos2piy, dy, cosh__2pi_xm1, sinh__2pi_xm1) +
        f_dx_l(dx + 2, cos2piy, dy, cosh__2pi_xp2, sinh__2pi_xp2) + f_dx_l(dx - 2, cos2piy, dy, cosh__2pi_xm2, sinh__2pi_xm2) +
        f_dx_l(dx + 3, cos2piy, dy, cosh__2pi_xp3, sinh__2pi_xp3) + f_dx_l(dx - 3, cos2piy, dy, cosh__2pi_xm3, sinh__2pi_xm3);
}

double nonIEEE_xy(double dx, double dy)
{
    double cosh2pix = cosh(M_PI * 2 * dx);
    double sinh2pix = sinh(M_PI * 2 * dx);
    double cos2piy = cos(M_PI * 2 * dy);
    if (dx < 0)
    {
        return
            f_l(dx + 5, cos2piy, dy, cosh__2pi_xp5) * (0 - dx) +
            f_l(dx - 4, cos2piy, dy, cosh__2pi_xm4) * (1 + dx) +
            f_l(dx + 4, cos2piy, dy, cosh__2pi_xp4) +
            g3(dx, cos2piy, dy, cosh2pix, sinh2pix);
    }

    return
        f_l(dx - 5, cos2piy, dy, cosh__2pi_xm5) * (0 + dx) +
        f_l(dx + 4, cos2piy, dy, cosh__2pi_xp4) * (1 - dx) +
        f_l(dx - 4, cos2piy, dy, cosh__2pi_xm4) +
        g3(dx, cos2piy, dy, cosh2pix, sinh2pix);
}

double nonIEEE_xy_diff_x(double dx, double dy)
{
    double cosh2pix = cosh(M_PI * 2 * dx);
    double sinh2pix = sinh(M_PI * 2 * dx);
    double cos2piy = cos(M_PI * 2 * dy);
    if (dx < 0)
    {
        return
            f_dx_l(dx + 5, cos2piy, dy, cosh__2pi_xp5, sinh__2pi_xp5) * (0 - dx) + f_l(dx + 5, cos2piy, dy, cosh__2pi_xp5) * (-1) +
            f_dx_l(dx - 4, cos2piy, dy, cosh__2pi_xm4, sinh__2pi_xm4) * (1 + dx) + f_l(dx - 4, cos2piy, dy, cosh__2pi_xm4) * (+1) +
            f_dx_l(dx + 4, cos2piy, dy, cosh__2pi_xp4, sinh__2pi_xp4) +
            g3_dx(dx, cos2piy, dy, cosh2pix, sinh2pix);
    }

    return
        f_dx_l(dx - 5, cos2piy, dy, cosh__2pi_xm5, sinh__2pi_xm5) * (0 + dx) + f_l(dx - 5, cos2piy, dy, sinh__2pi_xm5) * (+1) +
        f_dx_l(dx + 4, cos2piy, dy, cosh__2pi_xp4, sinh__2pi_xp4) * (1 - dx) + f_l(dx + 4, cos2piy, dy, sinh__2pi_xp4) * (-1) +
        f_dx_l(dx - 4, cos2piy, dy, cosh__2pi_xm4, sinh__2pi_xm4) +
        g3_dx(dx, cos2piy, dy, cosh2pix, sinh2pix);
}
#pragma endregion

int main()
{
#pragma region field calcualtion
    std::ofstream xy_difff("xy_difference.txt");
    std::ofstream IEEE_xyf("IEEE_xy.txt");
    std::ofstream nonIEEE_xyf("nonIEEE_xy.txt");
    std::ofstream f_f("f.txt");
    std::ofstream fl_f("fl.txt");
    std::ofstream cosh_f("cosh.txt");
    std::ofstream cosh_M_f("cosh_M.txt");

    std::ofstream gt_f("gt.txt");
    std::ofstream g3_f("g3.txt");

    if (!xy_difff || !g3_f || !gt_f || !f_f || !fl_f || !cosh_f || !cosh_M_f || !IEEE_xyf || !nonIEEE_xyf)
    {
        std::cerr << "Error: cannot open a txt file to write (field calculation).";
        exit(-1);
    }

    for (double y = 1. / 107183; y < 1; y += 1. / 100) // 107183 is prime, 1./107183 is almost 0. Numerically it is not good to use exactly 0.
    {
        {
            for (double x = 1. / 107183; x < 0.5; x += 1. / 100)
            {
                double dx = x;
                double dy = y;
                double cosh2pix = cosh(2 * M_PI * x);
                double sinh2pix = sinh(2 * M_PI * x);
                double cos2piy = cos(2 * M_PI * y);

                xy_difff << (IEEE_xy(x, y) - nonIEEE_xy(x, y)) / (IEEE_xy(x, y) + nonIEEE_xy(x, y)) << "\t";

                IEEE_xyf << IEEE_xy(x, y) << "\t";
                nonIEEE_xyf << nonIEEE_xy(x, y) << "\t";

                gt_f << g<3>(x, cos2piy, y) << "\t";
                g3_f << g3(x, cos2piy, y, cosh2pix, sinh2pix) << "\t";

                f_f << f(dx + ANALYTIC_FIELD_N, cos2piy, dy) << "\t";

                fl_f << f_l(dx + 4, cos2piy, dy, cosh__2pi_xp4) << "\t";

                cosh_f << cosh(2 * M_PI * (dx + 4)) << "\t";
                cosh_M_f << cosh2pix * cosh8pi + sinh2pix * sinh8pi << "\t";
            }

            xy_difff << "\n";
            IEEE_xyf << "\n";
            nonIEEE_xyf << "\n";

            gt_f << "\n";
            g3_f << "\n";

            f_f << "\n";
            fl_f << "\n";
            cosh_f << "\n";
            cosh_M_f << "\n";
        }
    }
#pragma endregion

#pragma region derivative field calculation
    std::ofstream xy_diff_xf("xy_diff_difference.txt");
    std::ofstream IEEE_xy_diff_xf("IEEE_xy_diff_x.txt");
    std::ofstream nonIEEE_xy_diff_xf("nonIEEE_xy_diff_x.txt");
    std::ofstream f_dx_f("f_dx.txt");
    std::ofstream g_dx_f("g_dx.txt");
    std::ofstream f_dx_Mf("f_dx_M.txt");
    std::ofstream g_dx_Mf("g_dx_M.txt");
    if (!xy_diff_xf || !IEEE_xy_diff_xf || !nonIEEE_xy_diff_xf || !f_dx_f || !g_dx_f || !f_dx_Mf || !g_dx_Mf)
    {
        std::cerr << "Error: cannot open a file to write (derivative calculation).";
        exit(-1);
    }

    for (double y = 1. / 107183; y < 1; y += 1. / 100) // 107183 is prime, 1./107183 is almost 0. Numerically it is not good to use exactly 0.
    {
        {
            for (double x = 1. / 107183; x < 1; x += 1. / 100)
            {
                double dx = x;
                double dy = y;
                double cos2piy = cos(2 * M_PI * y);
                double cosh2pix = cosh(2 * M_PI * x);
                double sinh2pix = sinh(2 * M_PI * x);

                xy_diff_xf << (IEEE_xy_diff_x(x, y) - nonIEEE_xy_diff_x(x, y)) / (IEEE_xy_diff_x(x, y) + nonIEEE_xy_diff_x(x, y)) << "\t";
                IEEE_xy_diff_xf << IEEE_xy_diff_x(x, y) << "\t";
                nonIEEE_xy_diff_xf << nonIEEE_xy_diff_x(x, y) << "\t";
                f_dx_f << f_dx(x - 5, cos2piy, y) << "\t";
                f_dx_Mf << f_dx_l(x - 5, cos2piy, y, cosh__2pi_xm5, sinh__2pi_xm5) << "\t";

                g_dx_f << g_dx<3>(x, cos2piy, y) << "\t";
                g_dx_Mf << g3_dx(x, cos2piy, y, cosh2pix, sinh2pix) << "\t";
            }

            xy_diff_xf << "\n";
            IEEE_xy_diff_xf << "\n";
            nonIEEE_xy_diff_xf << "\n";
            f_dx_f << "\n";
            g_dx_f << "\n";
            f_dx_Mf << "\n";
            g_dx_Mf << "\n";
        }
    }
#pragma endregion

    std::cout << "Done!\n";
}