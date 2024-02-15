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

#include "cosmology.h"
#include "survey.h"

// g++ -o pk2cgg -lgsl -lgslcblas pk2cgg.cc cosmology.cc survey.cc 

// Hubble reference parameter in km/s/Mpc
#define H0 100.

// speed of light in km/s
#define CSPEED 3.e+5

#define ZMAX 10.

using namespace std;

double cgg_integrand1(double x, void *p){

 struct f_pars * params = (struct f_pars *)p;

  double wL = (params->OmegaL);
  double wm = (params->Omegam);
  int band  = (params->band);
  int l     = (params->l);
  double h  = (params->h);
  
  // x ---> log10 of wavekumber k in units of Mpc^-1
  double k   = pow(10., x);
  double wg2 = pow(Wg(k, wL, wm, l, band, h), 2);
  
  return Delta2(k,h)*wg2;

}

double cgg_integrand2(double x, void *p){

 struct f_pars * params = (struct f_pars *)p;

  double wL = (params->OmegaL);
  double wm = (params->Omegam);
  int band  = (params->band);
  int l     = (params->l);
  double h  = (params->h);
  
  // x ---> comoving distance in units of Mpc

  double zval = r2z(x, wL, wm, h);
  double aval = 1./(1.+zval);
  double Da1 = growth(1., wL, wm);

  double kmin = KOH_MIN*h;
  double kmax = KOH_MAX*h;
  
  double ret;
  if ((l+0.5)/x>=kmin && (l+0.5)/x<=kmax)
    ret = pow(selection(band, zval)*growth(aval, wL, wm)/Da1/x,2)*PowerSpectrum((l+0.5)/x/h)*pow(Hubble(aval, wL, wm)*H0*h/CSPEED,2);
  else
    ret = 0.;
  
  return ret/h/h/h;

}

double cgg(double OmegaL, double Omegam, int l, int band, double h, double bg){

  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

  struct f_pars params;

  params.OmegaL = OmegaL;
  params.Omegam = Omegam;
  params.band = band;
  params.l = l;
  params.h = h;
  
  gsl_function F;

  double result, error;
  
  if (l<=20){
  
    F.function = &cgg_integrand1;

    F.params = &params;

    double logkmin = -4.+log10(h);
    double logkmax =  log10(2.)+log10(h);

    gsl_integration_qag (&F, logkmin, logkmax, 0., 1.e-5, 1000, 6, w, &result, &error);

    result *= 4.*M_PI*log(10.)*bg*bg;
    
  }
  else{

    F.function = &cgg_integrand2;

    F.params = &params;

    double rmin = 0.;
    double rmax =  z2r(ZMAX, OmegaL, Omegam, h);

    gsl_integration_qag (&F, rmin, rmax, 0., 1.e-5, 1000, 6, w, &result, &error);

    result *= bg*bg;

  }
    
  gsl_integration_workspace_free(w);
  
  return result;  

}


int main(int argc, char **argv){

  double OmegaL = 0.7; // cosmological constant to critical density ratio
  double Omegam = 0.3; // dark matter to critical density ratio
  double h = 0.7;      // Hubble parameter in units of 100 km/s/Mpc
  double bg = 1.; // galaxy bias
  
  string pkfname, outfname;
  int lmax = 128, band = 1;
  
  if (argc<15){
    cerr << "Usage: ./pk2cgg -pk [input power spectrum (ascii)] -wl [Omega Lambda] -wm [Omega matter] -h [H0 in km/s/Mpc] -lmax [maximum multipole] -band [selection function id] -bg [galaxy bias] -out [output C_gg file (ascii)]"<< endl;
    cerr << "Available selection functions:"<< endl;
    cerr << "-band 1 : Afshordi's parameterization of 2MASS (XSC) band 1 selection function."<< endl;
    cerr << "-band 2 : Afshordi's parameterization of 2MASS (XSC) band 2 selection function."<< endl;
    cerr << "-band 3 : Afshordi's parameterization of 2MASS (XSC) band 3 selection function."<< endl;
    cerr << "-band 4 : Afshordi's parameterization of 2MASS (XSC) band 4 selection function."<< endl;
    cerr << "-band 5 : Gaussian selection function."<< endl;
    cerr << "-band 6 : Spline to 2MASS (2MPZ) photo-z distribution for 12.0<K20<12.5."<< endl;
    cerr << "-band 7 : Spline to 2MASS (2MPZ) photo-z distribution for 12.5<K20<13.0."<< endl;
    cerr << "-band 8 : Spline to 2MASS (2MPZ) photo-z distribution for 13.0<K20<13.5."<< endl;
    cerr << "-band 9 : Spline to 2MASS (2MPZ) photo-z distribution for 13.5<K20<13.9."<< endl;
    cerr << "-band 10: Spline to 2MASS (2MPZ) photo-z distribution for 12.0<K20<13.9."<< endl;
    exit(1);
  }

  int idx;
  // read in user run flags and variables
  for (int i = 1; i < argc; i++) {
    idx = i;
    if (strstr(argv[i], "-pk" ) == argv[i] ) {i++; pkfname = argv[i];}   
    if (strstr(argv[i], "-out" ) == argv[i] ) {i++; outfname = argv[i];}
    if (strstr(argv[i], "-wl" ) == argv[i] ) {i++; OmegaL = atof(argv[i]);}
    if (strstr(argv[i], "-wm" ) == argv[i] ) {i++; Omegam = atof(argv[i]);}
    if (strstr(argv[i], "-bg" ) == argv[i] ) {i++; bg = atof(argv[i]);}
    if (strstr(argv[i], "-h" ) == argv[i] ) {i++; h = atof(argv[i]);}    
    if (strstr(argv[i], "-lmax" ) == argv[i] ) {i++; lmax = atoi(argv[i]);}    
    if (strstr(argv[i], "-band" ) == argv[i] ) {i++; band = atoi(argv[i]);} 
    
    if (idx==i){
      printf("Unknown option %s. Aborting. \n",argv[i]);
      exit(1);
    }
    
  }
  
  //string infname = "/Users/emoura/camb2015/camb/test_matterpower.dat";

  ofstream output(outfname.c_str());
  
  InitSpline(pkfname);

  if (band>=6 && band<=10)
    InitSelectionFunctionSpline(band);

  output << 0. << " " << 0. << endl; // monopole
  output << 0. << " " << 0. << endl; // dipole

  for (int l=2; l<=lmax; l++){

    double cggval = cgg(OmegaL, Omegam, l, band, h, bg);

    cerr << l << " " << cggval << endl;    
    output << l << " " << cggval << endl;
    
  }

  DestroySpline();

  if (band>=6 && band<=9)
    DestroySelectionSpline();
  
  output.close();
  
  return 0;

}
