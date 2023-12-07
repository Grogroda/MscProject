#include "cosmology.h"

double *x, *y;
gsl_interp_accel *acc;
gsl_spline *spline;


double *x2, *y2;
gsl_interp_accel *acc2;
gsl_spline *spline2;

// Hubble parameter in units of H0
double Hubble(double a, double OmegaL, double Omegam){

  return sqrt(OmegaL + Omegam/a/a/a + (1.-Omegam-OmegaL)/a/a);

}

double dHda(double a, double OmegaL, double Omegam){

  return 0.5*(-3.*Omegam/a/a/a/a - 2.*(1.-Omegam-OmegaL)/a/a/a)/Hubble(a, OmegaL, Omegam);

}
  
double growth_integrand(double x, void *p){
  
  struct f_pars * params = (struct f_pars *)p;

  // x ---> scale factor a
  double wL = (params->OmegaL);
  double wm = (params->Omegam);                                                                   

  return 1./pow(x*Hubble(x, wL, wm),3);

}
  
double growth(double a, double OmegaL, double Omegam){

  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  F.function = &growth_integrand;
  struct f_pars params;

  params.OmegaL = OmegaL;
  params.Omegam = Omegam;

  F.params = &params;

  double result, error;

  gsl_integration_qag (&F, 0., a, 0., 1.e-6, 1000, 6, w, &result, &error);

  gsl_integration_workspace_free(w);
  
  return 2.5*Omegam*Hubble(a, OmegaL,Omegam)*result;

}


double z2r_integrand(double x, void *p){
  
  struct f_pars * params = (struct f_pars *)p;

  double wL = params->OmegaL;
  double wm = params->Omegam;
  double h = params->h;
  
  double a = 1./(1.+x);
  
  return CSPEED/H0/h/Hubble(a, wL, wm);

}


double z2r(double z, double OmegaL, double Omegam, double h){

  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  F.function = &z2r_integrand;
  struct f_pars params;

  params.OmegaL = OmegaL;
  params.Omegam = Omegam;
  params.h = h;
  
  F.params = &params;

  double result, error;

  gsl_integration_qag (&F, 0., z, 0., 1.e-6, 1000, 6, w, &result, &error);

  gsl_integration_workspace_free(w);
  
  return result;

}

double z2r_for_root(double x, void *p){

  struct f_pars * params = (struct f_pars *)p;

  double wL   = params->OmegaL;
  double wm   = params->Omegam;
  double rval = params->rval;
  double h    = params->h;
  
  double ret = z2r(x, wL, wm, h) - rval;
  
  return ret;

}

double z2r_for_root_deriv(double x, void *p){

  struct f_pars * params = (struct f_pars *)p;
  
  double wL = params->OmegaL;
  double wm = params->Omegam;
  double h  = params->h;

  double ret = CSPEED/Hubble(1./(1.+x), wL, wm)/H0/h;    
  
  return ret;
  
}

void z2r_for_root_fdf(double x, void *p, double *y, double *dy){

  *y  = z2r_for_root(x, p); 
  *dy = z2r_for_root_deriv(x, p);

}


double r2z(double rval, double OmegaL, double Omegam, double h){

  if (rval==0.)
    return 0.;
  
  // need to implement
  bool verbo = false;
  
  int status;
  int iter = 0, max_iter = 100;

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  
  
  double x0, x = rval*ZL/z2r(ZL, OmegaL, Omegam, h);

  gsl_function F;
  
  struct f_pars params;
  params.OmegaL = OmegaL;
  params.Omegam = Omegam;
  params.rval   = rval;  
  params.h      = h;  
  
  F.function = &z2r_for_root;
  F.params = &params;
  
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, 0., 10.);

  if (verbo){
    printf ("%-5s %10s %10s\n",
            "iter", "root", "err(est)");

    printf ("%5s [%9s, %9s] %9s\n",
            "iter", "lower", "upper", "err(est)");
  }
  
  do
    {
      iter++;

      status = gsl_root_fsolver_iterate (s);
      x = gsl_root_fsolver_root (s);
      double x_lo = gsl_root_fsolver_x_lower (s);
      double x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,
                                       0, 0.001);

      if (verbo){
        if (status == GSL_SUCCESS)
          printf ("Converged:\n");
      }

      if (verbo){
        printf ("%5d [%.7f, %.7f] %.7f\n",
                iter, x_lo, x_hi, 
                x_hi - x_lo);
      }
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);
  
  return x;

}

double wg_integrand(double x, void *p){
  
  struct f_pars * params = (struct f_pars *)p;

  double wL = params->OmegaL;
  double wm = params->Omegam;
  double k  = params->k;
  double h  = params->h;
  double z0 = params->z0;
  double beta = params->beta;
  double lbda = params->lbda;
  int l     = params->l;

  // x ---> redshift z
  double a = 1./(1.+x);
  double cd = z2r(x, wL, wm, h);

  double Da1 = growth(1.,wL,wm);
  
  return selection(z0, beta, lbda, x)*gsl_sf_bessel_jl(l, k*cd)*growth(a,wL,wm)/Da1;

}


double wg_integrand_mc(double k, double z, int l, double wL, double wm, double h, double z0, double beta, double lbda){
  
  double a = 1./(1.+z);
  double cd = z2r(z, wL, wm, h);

  double Da1 = growth(1.,wL,wm);
  
  return selection(z0, beta, lbda, z)*gsl_sf_bessel_jl(l, k*cd)*growth(a,wL,wm)/Da1;

}

double Wg(double k, double OmegaL, double Omegam, int l, double z0, double beta, double lbda, double h){

  struct f_pars params;
  params.OmegaL = OmegaL;
  params.Omegam = Omegam;
  params.h      = h;
  params.z0     = z0;
  params.beta   = beta;
  params.lbda   = lbda;
  params.k      = k;
  params.l      = l;
  
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  F.function = &wg_integrand;

  F.params = &params;

  double zmax = ZMAX;
  double result, error;

  gsl_integration_qag (&F, 0., zmax, 0., 1.e-6, 1000, 6, w, &result, &error);

  gsl_integration_workspace_free(w);
  
  return result;

}

double wt_integrand(double x, void *p){
  
  struct f_pars * params = (struct f_pars *)p;

  double wL = params->OmegaL;
  double wm = params->Omegam;
  double k  = params->k;
  double h  = params->h;
  //  int band  = params->band;
  int l     = params->l;

  // x ---> ln of redshift z
  double z = exp(x);
  
  double a = 1./(1.+z);
  double cd = z2r(z, wL, wm, h);
  double Da1 = growth(1.,wL,wm);

  double f = pow(wm/pow(Hubble(a, wL, wm),2)/a/a/a, 0.6);

  double fac = -3.*wm*pow(H0*h/CSPEED,2)*TCMB/k/k;

  return fac*gsl_sf_bessel_jl(l, k*cd)*(f-1.)*(growth(a,wL,wm)/Da1)*z;

}


double wt_integrand_mc(double k, double lnz, int l, double wL, double wm, double h){

 double z = exp(lnz);
  
  double a = 1./(1.+z);
  double cd = z2r(z, wL, wm, h);
  double Da1 = growth(1.,wL,wm);

  double f = pow(wm/pow(Hubble(a, wL, wm),2)/a/a/a, 0.6);

  double fac = -3.*wm*pow(H0*h/CSPEED,2)*TCMB/k/k;

  return fac*gsl_sf_bessel_jl(l, k*cd)*(f-1.)*(growth(a,wL,wm)/Da1)*z;

}


double Wt(double k, double OmegaL, double Omegam, int l, double h){

  struct f_pars params;
  params.OmegaL = OmegaL;
  params.Omegam = Omegam;
  params.h      = h;
  params.k      = k;
  params.l      = l;

  
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  F.function = &wt_integrand;

  F.params = &params;

  double logzmin = log(0.0001);
  double logzmax = log(ZL);
  double result, error;

  gsl_integration_qag (&F, logzmin, logzmax, 0., 1.e-6, 1000, 6, w, &result, &error);

  gsl_integration_workspace_free(w);
  
  return result;

}



int ReadInputSpectrum(string fname){

  ifstream input(fname.c_str());

  vector<double> kvec, pvec;

  char line[2048];
  double koh, pk;

  input.getline(line, 2048);
  while (line[0]== '#') input.getline(line, 2048);

  while (!input.eof()){
 
    sscanf(line,"%lf %lf\n",&koh,&pk);
    kvec.push_back(koh);
    pvec.push_back(pk);
    input.getline(line,2048);
    if (input.eof())
      break;
  }

  int nks = pvec.size();
  // Comenting out this screen output
  // cerr << nks << " values of k read."<< endl;
  
  int i;
  x = (double *) malloc(nks*sizeof(double));
  y = (double *) malloc(nks*sizeof(double));

  //cerr << "log10(k [h Mpc^-1])    P(k) [h^-3 Mpc^3]" << endl;
  for (i = 0; i < nks; i++){
   x[i] = log10(kvec[i]);
   y[i] = pvec[i];
   // printf ("%15.6e       %15.6e\n", x[i], y[i]);
  }

  input.close();

  return nks;

}

void InitSpline(string fname){

  int nks = ReadInputSpectrum(fname);
  
  acc = gsl_interp_accel_alloc ();
  spline = gsl_spline_alloc (gsl_interp_cspline, nks);

  gsl_spline_init (spline, x, y, nks);
  
}

// units = h^-3 Mpc^3
double PowerSpectrum(double koh){

    //cerr << "PSepec "<< gsl_spline_eval (spline, log10(koh), acc) << endl;
  return gsl_spline_eval (spline, log10(koh), acc);
  
}


// adimensional
double Delta2(double k, double h){

  //cerr << "Delta2 " << PowerSpectrum(k/h) << endl;
  return 0.5*pow(k/h,3)*PowerSpectrum(k/h)/M_PI/M_PI;
  
}

void DestroySpline(){

  //cerr << "Spline destroyed" << endl;

  if (spline != NULL)
    gsl_spline_free (spline);

  if (acc != NULL)
    gsl_interp_accel_free (acc);

  if (x!=NULL)
    free(x);

  if(y!=NULL)
    free(y);
  
}


double WtSpline(int l, double k){

  return gsl_spline_eval (spline2, log10(k), acc2);
  
}

void InitWtSpline(int l, double OmegaL, double Omegam, double h){

  int npts = 500;

  if (x2==NULL)
    x2 = (double *) malloc(npts*sizeof(double));

  if (y2==NULL)
    y2 = (double *) malloc(npts*sizeof(double));

  double kmax = h*KOH_MAX;
  double kmin = h*KOH_MIN;
  double step = (log10(kmax)-log10(kmin))/(npts-1);
  
  //cerr << "InitWtSpline " << kmin << " "<< kmax << " " << step << endl;
  for (int i=0; i<npts; i++){
    x2[i] = log10(kmin)+i*step;
    y2[i] = Wt(pow(10., x2[i]), OmegaL, Omegam, l, h);
     cerr << i << " " << x2[i] << " " << y2[i] << endl;
  }
  
  acc2 = gsl_interp_accel_alloc ();
  spline2 = gsl_spline_alloc (gsl_interp_cspline, npts);

  gsl_spline_init (spline2, x2, y2, npts);
  
}

void DestroyWtSpline(){

  if (spline2 != NULL)
    gsl_spline_free (spline2);

  if (acc2 != NULL)
    gsl_interp_accel_free (acc2);

  if (x2!=NULL)
    free(x2);

  if(y2!=NULL)
    free(y2);
  
}

/*
double SphericalBessel(int l, double x){

  return gsl_spline_eval (spline2, log10(x), acc2);
  
}


void InitBesselSpline(int l, double OmegaL, double Omegam, double h){

  int npts = 200;

  if (x2==NULL)
    x2 = (double *) malloc(npts*sizeof(double));

  if (y2==NULL)
    y2 = (double *) malloc(npts*sizeof(double));



  mpfr_t x, y;
mpfr_init2 (x, 256);          // use default precision 
mpfr_init2 (y, 256);          // precision exactly 256 bits 

  double aux;
  long expo;
  
  double xmax = 2.* z2r(ZL, OmegaL, Omegam, h);
  double xmin = 1.e-6;
  double step = (log10(xmax)-log10(xmin))/npts;

  cerr << "InitBesselSpline " << xmin << " "<< xmax << " " << step << endl;
  for (int i=0; i<npts; i++){
    x2[i] = log10(xmin)+i*step;
    mpfr_init_set_d (x, pow(10.,x2[i]), GMP_RNDN);
    //    mpfr_jn (y, l, x, GMP_RNDN);
    mpfr_init_set_d (y, gsl_sf_bessel_jl(l, pow(10.,x2[i])), GMP_RNDN);
    y2[i] = mpfr_get_d (y, GMP_RNDN);
    //cerr << i << " " << x2[i] << " " << y2[i] << endl;
  }
  
  acc2 = gsl_interp_accel_alloc ();
  spline2 = gsl_spline_alloc (gsl_interp_cspline, npts);

  gsl_spline_init (spline2, x2, y2, npts);

// When the program is about to exit, do ... 
  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_free_cache ();
  
}


void DestroyBesselSpline(){

  if (spline2 != NULL)
    gsl_spline_free (spline);

  if (acc2 != NULL)
    gsl_interp_accel_free (acc);

  if (x2!=NULL)
    free(x);

  if(y2!=NULL)
    free(y);
  
}
*/
