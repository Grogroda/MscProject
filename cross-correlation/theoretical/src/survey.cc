#include "survey.h"

double *x_s, *y_s;
gsl_interp_accel *acc_s;
gsl_spline *spline_s;

double selection(double z0, double beta, double lbda, double z){

  double z0_gauss  = 1.0;
  double sig_gauss = 0.2;
  
  double window = 0.;

  if (z<=10){
    window  = beta*pow(z/z0,beta*lbda-1.);
    window *= exp(-pow(z/z0,beta));
    window /= gsl_sf_gamma(lbda);
    window /= z0;
  }
  else
    window = 0.;
    
  return window;
  
}

double selection2mpz(int band, double z){

  double ret = 0.;

  double zmax[5];

  zmax[0] = 0.1475;
  zmax[1] = 0.196;
  zmax[2] = 0.245;
  zmax[3] = 0.294;
  zmax[4] = 0.294;

  if (z<zmax[band-6]){
    //   cerr << z << " " << zmax[band-6] << endl;
    ret = gsl_spline_eval (spline_s, z, acc_s);
  }

  return ret;
  
}

int ReadRedshiftSelection(string fname){

  ifstream input(fname.c_str()); //arquivo de entrada do tipo ifstream com nome input

  vector<double> zvec, wvec;

  char line[2048]; //array de caracter, definindo a maior linha que pode ler
  double z, wg;
  while (!input.eof()){ //varre as linhas do arquivo
    input.getline(line,2048);
    if (input.eof()) //checa se chegou no final do arquivo
      break;
    sscanf(line,"%lf %lf\n",&z,&wg); //Quebra a linha, supondo que tem 2 variáveis do tipo double. %lf %lf\n (%ld para inteiro) indica que tem um real+espaço+real+quebra linha
    zvec.push_back(z); //Joga dentro de um vetor, memória alocada dinamicamente
    wvec.push_back(wg);
  }
  
  int nzs = wvec.size();
  cerr << nzs << " values of photo-z redshift bins read."<< endl;
  
  int i;
  x_s = (double *) malloc((nzs+1)*sizeof(double));
  y_s = (double *) malloc((nzs+1)*sizeof(double));

 
  cerr << "photo-z       Selection function" << endl;
  for (i = 0; i < nzs+1; i++){
    if (i==0){
      x_s[0] = 0.;
      y_s[0] = 0.; 
    }
    else{
      x_s[i] = zvec[i];
      y_s[i] = wvec[i];
    }
    printf ("%15.6e       %15.6e\n", x_s[i], y_s[i]);
  }

  input.close();

  return nzs;

}

void InitSelectionFunctionSpline(int band){

  int nzs;

  //  string fname = ;
  
  if (band == 6)
    nzs = ReadRedshiftSelection("/Users/emoura/physics/cosmos/data/2mass_2mpz/photoz_selection_2mpz_12p0_k20_12p5.dat");
  else if (band == 7)
    nzs = ReadRedshiftSelection("/Users/emoura/physics/cosmos/data/2mass_2mpz/photoz_selection_2mpz_12p5_k20_13p0.dat");
  else if (band == 8)
    nzs = ReadRedshiftSelection("/Users/emoura/physics/cosmos/data/2mass_2mpz/photoz_selection_2mpz_13p0_k20_13p5.dat");
  else if (band == 9)
    nzs = ReadRedshiftSelection("/Users/emoura/physics/cosmos/data/2mass_2mpz/photoz_selection_2mpz_13p5_k20_13p9.dat");
  else if (band == 10)
    nzs = ReadRedshiftSelection("/Users/emoura/physics/cosmos/data/2mass_2mpz/photoz_selection_2mpz_12p0_k20_13p9.dat");
  else{
    cerr << "Unknown 2MPZ band. Aborting!"<< endl;
    exit(1);
  }
    
  
  acc_s = gsl_interp_accel_alloc ();
  spline_s = gsl_spline_alloc (gsl_interp_cspline, nzs);

  gsl_spline_init (spline_s, x_s, y_s, nzs);
  
}


void DestroySelectionSpline(){

  if (spline_s != NULL)
    gsl_spline_free (spline_s);

  if (acc_s != NULL)
    gsl_interp_accel_free (acc_s);

  if (x_s!=NULL)
    free(x_s);

  if(y_s!=NULL)
    free(y_s);
  
}

