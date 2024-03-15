#ifndef _COSMOS_H_
#define _COSMOS_H_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

//#include <gmp.h>
//#include <mpfr.h>

#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"

#include "survey.h"

// Hubble reference parameter in km/s/Mpc
#define H0 100.

// speed of light in km/s
#define CSPEED 3.e+5

#define ZMIN 1e-4
#define ZMAX 10.
#define ZL 1089.

#define KOH_MIN 1.e-4
#define KOH_MAX 2.

// CMB temperature today in microK
#define TCMB 2.726e+6


using namespace std;


struct f_pars {
  double OmegaL, Omegam, k, h, rval, z0, beta, lbda;
  int l;
};

//Beginning
extern "C"{
double Hubble(double a, double OmegaL, double Omegam);

double dHda(double a, double OmegaL, double Omegam);
  
double growth_integrand(double x, void *p);

double growth(double a, double OmegaL, double Omegam);

double z2r_integrand(double x, void *p);

double z2r(double z, double OmegaL, double Omegam, double h);

double z2r_for_root(double x, void *p);

double z2r_for_root_deriv(double x, void *p);

void z2r_for_root_fdf(double x, void *p, double *y, double *dy);

double r2z(double rval, double OmegaL, double Omegam, double h);

double wg_integrand(double x, void *p);

double wg_integrand_mc(double k, double z, int l, double wL, double wm, double h, double z0, double beta, double lbda);

double Wg(double k, double OmegaL, double Omegam, int l, double z0, double beta, double lbda, double h);

double wt_integrand(double x, void *p);

double wt_integrand_mc(double k, double z, int l, double wL, double wm, double h);

double Wt(double k, double OmegaL, double Omegam, int l, double h);

int DefineSpectrum(double karr[], double pkarr[], int nks);

int ReadInputSpectrum(string fname);

void InitSpline(bool read_input, string fname, double karr[], double pkarr[], int nks);

// units = h^-3 Mpc^3
double PowerSpectrum(double koh);

// adimensional
double Delta2(double k, double h);

void DestroySpline();
/*
double SphericalBessel(int l, double x);
void InitBesselSpline(int l, double OmegaL, double Omegam, double h);
void DestroyBesselSpline();
*/

double WtSpline(int l, double x);
void InitWtSpline(int l, double OmegaL, double Omegam, double h);
void DestroyWtSpline();
}
//End

#endif
