#include <iostream>
#include <fstream>
#include <string>

#include "correlations.h"

#include <cmath>
#include <math.h> 
#include <gsl/gsl_multimin.h>

#define LBREAK 10

#define _USE_MATH_DEFINES

//Passar as funções todas para um novo código sem main()

/*  
Raciocinando:

Tenho que receber (z0,beta,lbda) -> Determinar Cltg para esses valores (0 <= l <= 128)
-> Calcular logP=sum(log(gauss(cltg))) -> Retornar P=exp(logP)
*/

using namespace std;

struct gaussians{
	vector<double> muvec;
	vector<double> varvec;
};

struct min_pars{
	double Omegam; //Mass energy density relative to critical mass
	double OmegaL; //Lambda energy density relative to critical mass
	double h; //Hubble parameter divided by 100
	double bg; //Bias factor for the selection function
	
	int mode; //mode==0 -> Gaussian quadrature integration; else -> MC integration
	int ncalls; //Number of MC calls
	
	string pkname; //Name of 3d matter power spectrum file name
	string gauss_name; //Name of gaussians table file name
	
	vector<double> ctg2MASS; //Cross-correlation power spectrum vector for a 2MASS band
};

double gauss(double x, double mu, double var){
	/*
	This function simply calculates a gaussian in the point x
	for given average and variance.
	*/

	double N = 1/(var*sqrt(2*M_PI));
	double expo = exp(pow((x-mu),2)/(2*var));

	return N*expo;	
}

gaussians read_gaussians(string fname){

	/*
	This function reads the gaussians file and returns
	an instance of the gaussian structure.
	*/

	ifstream input(fname.c_str()); //Input file of the type ifstream named input
	
	gaussians G;
	
	vector<int> lvec;
	vector<double> muvec, varvec;

	char line[2048]; //Array of characters, with an arbitrary limit to the number of lines
	
	int l;
	double mu, var;
	
	input.getline(line,2048);
	while (!input.eof()){ //Looks through every line
		
		sscanf(line,"%ld %lf %lf\n",&l,&mu,&var); //Breaks the line, assuming 3 columns: %ld %lf %lf\n=int+ +double+ +double+linebreak
	  	lvec.push_back(l); //Passes all values to their respective vectors
	  	muvec.push_back(mu);
	  	varvec.push_back(var);
	  	
	  	//cerr << "mu=" << mu << " & var=" << var << endl;
	  	
		input.getline(line,2048);
		if (input.eof()) //checa se chegou no final do arquivo
			break;
	}
	
	//Sets the structure's vectors:
	G.muvec = muvec;
	G.varvec = varvec;
	
	return G;
}

vector<double> read_2MASS(string fname){
	
	ifstream input(fname.c_str()); //arquivo de entrada do tipo ifstream com nome input
	
	vector<double> ctgvec;

	char line[2048]; //array de caracter, definindo a maior linha que pode ler

	double cltg;
	int l;
	
	input.getline(line,2048);
	while (!input.eof()){ //varre as linhas do arquivo
		
		sscanf(line,"%ld %lf\n",&l,&cltg); //Quebra a linha, supondo que tem 2 variáveis do tipo double. %lf %lf\n (%ld para inteiro) indica que tem um real+espaço+real+quebra linha
	  	ctgvec.push_back(cltg); //Joga dentro de um vetor, memória alocada dinamicamente
	  	
		input.getline(line,2048);
		if (input.eof()) //checa se chegou no final do arquivo
			break;
	}
	
	return ctgvec;
}

double calc_logP(vector<double> Ctg, gaussians G){
	/*Given a cross-correlation power spectrum and an
	instance of a gaussians structure, calculates and
	returns the log(P)=sum(log(p)).
	*/
	
	vector<double> muvec = G.muvec;
	vector<double> varvec = G.varvec;

	double logp = 0.0;
	for (int l = 2; l < 129; l++){
		logp+=log(gauss(Ctg.at(l), muvec.at(l), varvec.at(l))); //tirar o log
	}

	return logp;
} 


double P_relative(double beta, double z0, double lbda, double OmegaL, double Omegam, double h, double bg, int mode, int ncalls, string pkfname, gaussians Gs, vector<double> ctg2mass){
	
	InitSpline(pkfname);
	
	vector<double> CTG;
	
	//Monopole and dipole terms will not be used
	CTG.push_back(0.0); //l=0
	CTG.push_back(0.0); //l=1
	
	//#pragma omp parallel for
	for (int l=2; l<=LBREAK; l++){
		cerr << "l=" << l << endl;
		
		//#pragma omp critical
		//{
		double ctgval = ctg(OmegaL, Omegam, l, z0, beta, lbda, h, bg, mode, ncalls);
		
		CTG.push_back(ctgval);
		//}

		cerr << l << " " << ctgval << endl; // write cross-correlation to screen
		//output << l << " " << ctgval << endl; // write cross-correlation to output file
	}


	// Loop over the remaining multipoles (l>10)
	for (int l=LBREAK+1; l<129; l++){

		double ctgval = ctg(OmegaL, Omegam, l, z0, beta, lbda, h, bg, mode, ncalls);
		CTG.push_back(ctgval);

		cerr << l << " " << ctgval << endl; // write cross-correlation to screen
		//output << l << " " << ctgval << endl; // write cross-correlation to output file
	}
	
	double logp = calc_logP(CTG, Gs);
	cout << "log(P)=" << logp << endl;
	
	double logp2mass = calc_logP(ctg2mass, Gs);
	
	cout << "log(P)-log(P2mass)=" << logp-logp2mass << endl;
	
	if(beta < 0 || lbda < 0 || z0 < 0){
		return 1000*exp(logp-logp2mass);
	} else{
		return exp(logp-logp2mass);
	}
}



double Pgsl(const gsl_vector *v, void *params){

	double beta, z0, lbda;
	min_pars *p = (min_pars *)params;
	
	double OmegaL           = p->OmegaL;
	double Omegam           = p->Omegam;
	double h                = p->h;
	double bg               = p->bg;
	
	int mode                = p->mode;
	int ncalls              = p->ncalls;
	
	string pk               = p->pkname;
	string gauss_fname      = p->gauss_name;
	
	vector<double> ctg2mass = p->ctg2MASS;

	beta = gsl_vector_get(v, 0);
	z0 = gsl_vector_get(v, 1);
	lbda = gsl_vector_get(v, 2);
	
	gaussians Gs = read_gaussians(gauss_fname); 
	
	double Prob = P_relative(beta, z0, lbda, OmegaL, Omegam, h, bg, mode, ncalls, pk, Gs, ctg2mass);
	
	return Prob;
} 

//double *p = (double *)params; -> casting, converte do tipo void para um tipo determinado pelo usuario

/*
int main(){
  
  double OmegaL = 0.7;
  double Omegam = 0.3;
  double h = 0.7;
  double bg = 1.0;
  int mode = 1;
  int ncalls = 1000000;
  
  double beta = 1.0;
  double z0 = 0.05;
  double lbda = 5.0;
  
  string pk = "Cl_Pk/pk_3dmatter.dat";
  
  cerr << "Começando o cálculo:" << endl;
  
  double logprob=logP(beta, z0, lbda, OmegaL, Omegam, h, bg, mode, ncalls, pk);
  
  cout << "logP(1.0,0.05,5.0)=" << logprob << endl;
  
  return 0;
} 
*/


int main(void){
  
  min_pars pars;
  
  pars.OmegaL     = 0.7;
  pars.Omegam     = 0.3;
  pars.h          = 0.7;
  pars.bg         = 1.0;
  pars.mode       = 1;
  pars.ncalls     = 1000000;
  pars.pkname     = "tables/pk_3dmatter.dat"; //Generalizar depois
  pars.gauss_name = "tables/hists/gauss_fits.txt"; //Generalizar depois
  
  string fname2mass = "tables/2MASS/band1_bg1.37.dat"; //Generalizar depois
  
  //The paths here are relative to the main directory cross-correlation
  //I'll work on providing a user option in the future
  
  vector<double> ctg2massvec = read_2MASS(fname2mass);
  
  pars.ctg2MASS = ctg2massvec;

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  // Starting point

  x = gsl_vector_alloc(3);
  gsl_vector_set(x, 0, 1.0);
  gsl_vector_set(x, 1, 0.05);
  gsl_vector_set(x, 2, 5.0);
  
  // Set initial step sizes to 1 

  ss = gsl_vector_alloc(3);
  gsl_vector_set(ss, 0, 0.1);
  gsl_vector_set(ss, 1, 0.01);
  gsl_vector_set(ss, 2, 0.2);

  // Initialize method and iterate 
  minex_func.n = 3;
  minex_func.f = Pgsl; 
  minex_func.params = &pars;

  s = gsl_multimin_fminimizer_alloc (T, 3);
  
  cerr << "Memory Allocated. Setting minimizer next." << endl;
  
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss); 
  /*
  Apparently, for the algorithm I'm using, the above function
  will call my minimizing function N+1=4 times to setup the
  minimizer before it actually starts iterating, and since it's
  a slow function to calculate each point, it can take a while 
  for the program to start the actual iteration.  
  */
  
  cerr << "About to start loop." << endl; 

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      cerr << "Inside the loop!!" << endl;

      if (status)
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-2);

      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }

      printf ("%5d %10.3e %10.3e %10.3e f() = %7.3f size = %.3f\n",
              iter,
              gsl_vector_get (s->x, 0),
              gsl_vector_get (s->x, 1),
              gsl_vector_get (s->x, 2),
              s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < 100);

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return status;
} 
