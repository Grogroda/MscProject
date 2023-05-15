#ifndef _SURVEY_H_
#define _SURVEY_H_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"

using namespace std;

double selection(double z0, double beta, double lbda, double z);
int ReadRedshiftSelection(string fname);
void InitSelectionFunctionSpline(int band);
void DestroySelectionSpline();
double selection2mpz(int band, double z);

#endif
