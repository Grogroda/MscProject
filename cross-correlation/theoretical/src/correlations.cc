// Hubble reference parameter in km/s/Mpc (the actual value is a fraction of this value)
#define H0 100.

// speed of light in km/s
#define CSPEED 3.e+5

#define ZMAX 10.

#include "correlations.h"
#include "survey.h"
#include <string>

using namespace std;


void display_results (char *title, double result, double error)
//Comenting out the screen output
{
  //  printf ("%s ==================\n", title);
  //printf ("result = % .6f\n", result);
  //printf ("sigma  = % .6f\n", error);
}

double ctt_integrand(double x, void *p){
	struct f_pars * params = (struct f_pars *)p;

	double wL = (params->OmegaL);
	double wm = (params->Omegam);
	int l     = (params->l);
	double h  = (params->h);

	/*
	double z0   = (params->z0);
	double beta = (params->beta);
	double lbda = (params->lbda);
	*/

	//x ---> log10 of wavekumber k in units of Mpc^-1
	
	double k  = pow(10.,x);
	double wt = WtSpline(l, k);  
	double wt2 = pow(wt, 2);

	return Delta2(k,h)*wt2;
}


double cgg_integrand1(double x, void *p){
	struct f_pars * params = (struct f_pars *)p;

	double wL = (params->OmegaL);
	double wm = (params->Omegam);
	int l     = (params->l);
	double h  = (params->h);

	double z0   = (params->z0);
	double beta = (params->beta);
	double lbda = (params->lbda);

	//x ---> log10 of wavekumber k in units of Mpc^-1
	
	double k  = pow(10.,x);
	double wg2 = pow(Wg(k, wL, wm, l, z0, beta, lbda, h), 2);

	return Delta2(k,h)*wg2;

}

double cgg_integrand2(double x, void *p){

	struct f_pars * params = (struct f_pars *)p;

	double wL = (params->OmegaL);
	double wm = (params->Omegam);
	int l     = (params->l);
	double h  = (params->h);
	
	double z0   = (params->z0);
	double beta = (params->beta);
	double lbda = (params->lbda);

	// x ---> comoving distance in units of Mpc
	
	double zval = r2z(x, wL, wm, h);
	double aval = 1./(1.+zval);
	double Da1 = growth(1., wL, wm);

	double kmin = KOH_MIN*h;
	double kmax = KOH_MAX*h;

	double ret;

	if ((l+0.5)/x>=kmin && (l+0.5)/x<=kmax)
		ret = pow(selection(z0, beta, lbda, zval)*growth(aval, wL, wm)/Da1/x,2)*PowerSpectrum((l+0.5)/x/h)*pow(Hubble(aval, wL, wm)*H0*h/CSPEED,2);
	else
		ret = 0.;

	return ret/h/h/h;

}


double ctg_integrand1(double x, void *p){

 struct f_pars * params = (struct f_pars *)p;

  double wL = params->OmegaL;
  double wm = params->Omegam;
  double z0 = params->z0;
  double beta = params->beta;
  double lbda = params->lbda;
  int l     = params->l;
  double h  = params->h;
  
  // x ---> log10 of wavekumber k in units of Mpc^-1
  double k   = pow(10., x);
  double wg = Wg(k, wL, wm, l, z0, beta, lbda, h);

  //  double wt = Wt(k, wL, wm, l, h);
  double wt = WtSpline(l, k);  
  
  return Delta2(k,h)*wt*wg;

}


double ctg_integrand2(double x, void *p){

 struct f_pars * params = (struct f_pars *)p;

  double wL = params->OmegaL;
  double wm = params->Omegam;
  double z0 = params->z0;
  double beta = params->beta;
  double lbda = params->lbda;
  int l     = params->l;
  double h  = params->h;
  
  double zval = x;

  // turn redshift into comoving distance
  double r = z2r(zval, wL, wm, h);
  
  double aval = 1./(1.+zval);
  double Da1 = growth(1., wL, wm);

  // Attention: this is valid only for LCDM
  double f = pow(wm/pow(Hubble(aval, wL, wm),2)/aval/aval/aval, 0.6);
  
  if ((l+0.5)/r/h<1.e-4 || (l+0.5)/r/h>2.)
    return 0.;

  double ret = selection(z0, beta, lbda, zval)*pow(growth(aval, wL, wm)/Da1,2)*Hubble(aval, wL, wm)*(f-1.)*PowerSpectrum((l+0.5)/r/h);
  
  
  return ret/h/h/h;

}


double ctg_integrand3(double *x, size_t dim, void *p){

 struct f_pars * params = (struct f_pars *)p;

  double wL = params->OmegaL;
  double wm = params->Omegam;
  double z0 = params->z0;
  double beta = params->beta;
  double lbda = params->lbda;
  int l     = params->l;
  double h  = params->h;
  
  // x ---> log10 of wavekumber k in units of Mpc^-1
  double k    = pow(10., x[0]);
  double z    = x[1];
  double lnzp = x[2];  
  //double wg = Wg(k, wL, wm, l, z0, beta, lbda, h);

  //  double wt = Wt(k, wL, wm, l, h);
  //double wt = WtSpline(l, k);  

  //cerr << "integ3" << " " << Delta2(k,h) << " " << wt_integrand_mc(k, lnzp, l, wL, wm, h) << " " << wg_integrand_mc(k, z, l, wL, wm, h, z0, beta, lbda) << endl;

  return Delta2(k,h)*wt_integrand_mc(k, lnzp, l, wL, wm, h)*wg_integrand_mc(k, z, l, wL, wm, h, z0, beta, lbda);

}

double ctt(double OmegaL, double Omegam, int l, double h, double bg){

  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
  
  struct f_pars params;

  params.OmegaL = OmegaL;
  params.Omegam = Omegam;
  /*
  params.z0 = z0;
  params.beta = beta;
  params.lbda = lbda;
  */
  params.l = l;
  params.h = h;

  double result, error;

  gsl_function F;

  F.function = &ctt_integrand;
  F.params = &params;

  double logkmin = log10(KOH_MIN*h);
  double logkmax = log10(KOH_MAX*h);

  //cerr << "Starting integration" << endl;

  gsl_integration_qag (&F, logkmin, logkmax, 0., 1.e-6, 1000, 6, w, &result, &error);

  //cerr << "Integration finished" << endl;

  result *= 4.*M_PI*log(10.)*bg;
  
  gsl_integration_workspace_free(w);

  return result;
  
}


double cgg(double OmegaL, double Omegam, int l, double z0, double beta, double lbda, double h, double bg){
	
	gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

	struct f_pars params;

	params.OmegaL = OmegaL;
	params.Omegam = Omegam;
	params.l      = l;
	params.h      = h;

	params.z0   = z0;
	params.beta = beta;
	params.lbda = lbda;

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
		double rmax = z2r(ZMAX, OmegaL, Omegam, h);

		gsl_integration_qag (&F, rmin, rmax, 0., 1.e-5, 1000, 6, w, &result, &error);

		result *= bg*bg;

	}

	gsl_integration_workspace_free(w);

	return result;

}


double ctg_quad(double OmegaL, double Omegam, int l, double z0, double beta, double lbda, double h, double bg){

  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

  struct f_pars params;

  params.OmegaL = OmegaL;
  params.Omegam = Omegam;
  params.z0 = z0;
  params.beta = beta;
  params.lbda = lbda;
  params.l = l;
  params.h = h;

  double result, error;
  
  gsl_function F;

  if (l<=10){
  
    F.function = &ctg_integrand1;
    F.params = &params;

    double logkmin = log10(KOH_MIN*h);
    double logkmax = log10(KOH_MAX*h);

    gsl_integration_qag (&F, logkmin, logkmax, 0., 1.e-6, 1000, 6, w, &result, &error);

    result *= 4.*M_PI*log(10.)*bg;

  }
  else{ // if l>=10, we will use the Limber approximation for the spherical Bessel functions

    F.function = &ctg_integrand2;
    F.params = &params;
    
    double zmin = 0.0;
    double zmax = ZMAX;
    gsl_integration_qag (&F, zmin, zmax, 0., 1.e-6, 1000, 6, w, &result, &error);
    
    result *= -3*bg*TCMB*Omegam*pow(H0*h/CSPEED,3)/pow(l+0.5,2);
    
  }

  gsl_integration_workspace_free(w);  
  
  return result;  

}


double ctg_mc(double OmegaL, double Omegam, int l, double z0, double beta, double lbda, double h, double bg, int ncalls){

  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

  struct f_pars params;

  params.OmegaL = OmegaL;
  params.Omegam = Omegam;
  params.z0 = z0;
  params.beta = beta;
  params.lbda = lbda;
  params.l = l;
  params.h = h;

  double result, error;

  double logkmin = log10(KOH_MIN*h);
  double logkmax = log10(KOH_MAX*h);
  
  double xl[3] = { logkmin, ZMIN, log(ZMIN)};
  double xu[3] = { logkmax, ZMAX, log(ZL)  };

  gsl_function F;

  const gsl_rng_type *T;
  gsl_rng *r;
  
  gsl_monte_function G = { &ctg_integrand3, 3, &params };

  size_t calls = ncalls;

  if (l<=10){

    gsl_rng_env_setup ();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    {
      gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

      gsl_monte_vegas_integrate (&G, xl, xu, 3, 1000, r, s, &result, &error);

      display_results ("vegas warm-up", result, error);

      //Comenting out screen output

      //      printf ("converging...\n");

      do
	{
	  gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s, &result, &error);
	  //printf ("result = % .6f sigma = % .6f "
	  //	  "chisq/dof = %.1f\n", result, error, gsl_monte_vegas_chisq (s));
	}
      while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
      
      display_results ("vegas final", result, error);

      gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);
   

    result *= 4.*M_PI*log(10.)*bg;

  }
  else{ // if l>=10, we will use the Limber approximation for the spherical Bessel functions

    F.function = &ctg_integrand2;
    F.params = &params;
    
    double zmin = 0.0;
    double zmax = ZMAX;
    gsl_integration_qag (&F, zmin, zmax, 0., 1.e-6, 1000, 6, w, &result, &error);
    
    result *= -3*bg*TCMB*Omegam*pow(H0*h/CSPEED,3)/pow(l+0.5,2);
    
  }

  gsl_integration_workspace_free(w);  
  
  return result;  

}

double cgg4py(double OmegaL, double Omegam, int l, double z0, double beta, double lbda, double h, double bg, int mode, int ncalls, double karr[], double pkarr[]){

	InitSpline(read_input = false, karr = karr, pkarr = pkarr);

	double cggl=cgg(OmegaL, Omegam, l, z0, beta, lbda, h, bg);
	return cggl;
}

double ctg4py(double OmegaL, double Omegam, int l, double z0, double beta, double lbda, double h, double bg, int mode, int ncalls, double karr[], double pkarr[]){
	
	InitSpline(read_input = false, karr = karr, pkarr = pkarr)

	double ctgl=ctg(OmegaL, Omegam, l, z0, beta, lbda, h, bg, mode, ncalls);
	return ctgl;
}


double ctg(double OmegaL, double Omegam, int l, double z0, double beta, double lbda, double h, double bg, int mode, int ncalls){

  double ret;
  
  if (mode == 0){
    ret = ctg_quad(OmegaL, Omegam, l, z0, beta, lbda, h, bg);
    }
  else{
    ret = ctg_mc(OmegaL, Omegam, l, z0, beta, lbda, h, bg, ncalls);
    }

  return ret;

}
