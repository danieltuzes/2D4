
// PeriodicShearStressELTE.h : contains function declarations for PeriodicShearStressELTE.cpp, also included in project_parser.cpp

/*
 * NOTE: This class is based on a library which was maintained by the Department of Materials physics
 * at Eötvös Loránd University, Budapest
 */

#ifndef SDDDST_CORE_PERIODIC_SHEAR_STRESS_H
#define SDDDST_CORE_PERIODIC_SHEAR_STRESS_H

#define VERSION_periodic_shear_stress_elte 0.1

#include <vector>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>

//#include <gsl/gsl_sf_bessel.h>

#include "Fields/Field.h"

/*#define STRESSMSIZE	1024*/
#define STRESSMDT	double
//#define STRESSMFNXX	"/home/ispanovity/work/progs/dislocation_dynamics/periodic_sigma/periodic_stress_xx_1024x1024_bin.dat"
//#define STRESSMFNXX	"periodic_stress_xx_1024x1024_bin.dat"
//#define STRESSMFNYY	"/home/ispanovity/work/progs/dislocation_dynamics/periodic_sigma/periodic_stress_yy_1024x1024_bin.dat"
//#define STRESSMFNYY	"periodic_stress_yy_1024x1024_bin.dat"
//#define STRESSMFNXY	"/home/ispanovity/work/progs/dislocation_dynamics/periodic_sigma/periodic_stress_xy_1024x1024_bin.dat"
//#define STRESSMFNXY	"periodic_stress_xy_1024x1024_bin.dat"
//#define STRESSMFNXYCUT	"/home/ispanovity/work/progs/dislocation_dynamics/periodic_sigma/periodic_stress_xy_cut_1024x1024_bin.dat"
//#define STRESSMFNXY2ND	"/home/ispanovity/work/progs/dislocation_dynamics/periodic_sigma/periodic_stress_xy_2nd_order_1024x1024_bin.dat"
/*#define STRESSMRADIN2	(1000.0/STRESSMSIZE/STRESSMSIZE)*/
/*#define STRESSMRADOUT2	(1.3*STRESSMRADIN2)*/
#define STRESSMFACT	(4 * M_PI)

#define MAX_H_FOR_DIFF (1.0/(4*stress_matrix_size))
#define MAX_DELTA_FOR_DIFF 1e-10
#define MIN_H_FOR_DIFF (MAX_H_FOR_DIFF/1000000)

#define ERRNO_INPUT	1

#define ANHARM_R0 0.001

namespace sdddstCoreELTE {

class PeriodicShearStressELTE : public sdddstCore::Field
{
public:	
    PeriodicShearStressELTE();
    virtual ~PeriodicShearStressELTE();
	
    void loadStress(std::string path, const char* str, int n, double R0 = 0);

    virtual double xy(double x, double y);
    virtual double xy_diff_x(double x, double y);

    bool get_field_loaded() const;

    bool get_diff_x_field_loaded() const;

protected:
    double** stressm_xy;
    double** stressm_xy_diff_x;

    double r0;
    unsigned int stress_matrix_size;
    double stress_interp_radius_in_sq;
    double stress_interp_radius_out_sq;

    bool field_loaded;
    bool diff_x_field_loaded;
};

}
#endif //PERIODIC_SHEAR_STRESS_H
