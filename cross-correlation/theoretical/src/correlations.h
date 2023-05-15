#ifndef _CORRELATIONS_H_
#define _CORRELATIONS_H_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>

#include "cosmology.h"
#include <cstring>
#include "survey.h"

using namespace std;

void display_results (char *title, double result, double error);

double ctg_integrand1(double x, void *p);

double ctg_integrand2(double x, void *p);

double ctg_integrand3(double *x, size_t dim, void *p);

double ctg_quad(double OmegaL, double Omegam, int l, double z0, double beta, double lbda, double h, double bg);

double ctg_mc(double OmegaL, double Omegam, int l, double z0, double beta, double lbda, double h, double bg, int ncalls);

double ctg(double OmegaL, double Omegam, int l, double z0, double beta, double lbda, double h, double bg, int mode, int ncalls);

#endif
