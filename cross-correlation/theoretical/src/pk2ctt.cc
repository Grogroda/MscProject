#include "correlations.h"
#include <omp.h>
#include <chrono>
#include <thread>

using namespace std;

// g++ -o pk2ctg -lgsl -lgslcblas pk2ctg.cc cosmology.cc survey.cc 

#define LBREAK 10

// Given an input 3D matter power spectrum for a LCDM model, this code calculates the expected angular cross-correlation power spectrum
// in harmonic space between the CMB temperature and the galaxy contrast maps.

//Function to display time for testing:
string getTimeStr(){
    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    std::string s(30, '\0');
    std::strftime(&s[0], s.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    return s;
}

int main(int argc, char **argv){

  // default valeus, can all be changed at command line
  double OmegaL = 0.7; // cosmological constant to critical energy density ratio 
  double Omegam = 0.3; // dark matter to critical energy density ratio
  double h = 0.7;      // Hubble parameter in units of 100 km/s/Mpc
  double bg = 1.;      // galaxy bias
  /*
  double z0 = 0.054;    // survey parameter z0
  double beta = 1.8;    // survey parameter beta
  double lbda = 1.6;    // survey parameter lambda
  */
  
  string pkfname, outfname;
  int lmax = 128, band = 1;
  
  if (argc<6){
    cerr << "Usage: ./pk2ctg -pk [input power spectrum (ascii)] -wl [Omega Lambda] -wm [Omega matter] -h [H0 in km/s/Mpc] -lmax [maximum multipole] -bg [galaxy bias] -out [output C_tg file (ascii)]" << endl;
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

  /*
  int mode = 0; // 0 -> quadrature, != 0 -> Monte Carlo
  int ncalls = 500000;
  int nthreads = -1;
  */
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
    /*
    if (strstr(argv[i], "-z0" ) == argv[i] ) {i++; z0 = atof(argv[i]);} 
    if (strstr(argv[i], "-beta" ) == argv[i] ) {i++; beta = atof(argv[i]);} 
    if (strstr(argv[i], "-lbda" ) == argv[i] ) {i++; lbda = atof(argv[i]);}
    if (strstr(argv[i], "-mode" ) == argv[i] ) {i++; mode = atoi(argv[i]);}    
    if (strstr(argv[i], "-ncalls" ) == argv[i] ) {i++; ncalls = atoi(argv[i]);}
    if (strstr(argv[i], "-nthreads" ) == argv[i] ) {i++; nthreads = atoi(argv[i]);} 
    */
    
    if (idx==i){
      printf("Unknown option %s. Aborting. \n",argv[i]);
      exit(1);
    }
    
  }

/*
  if (nthreads > 0){
    omp_set_dynamic(0);
    omp_set_num_threads(nthreads);
    cerr << nthreads << " thread(s) will be used." << endl;
  }
*/

  cerr << "Matter 3D power spectrum will be read from " << pkfname << endl;
  cerr << "Output file will be written to " << outfname << endl;
  cerr << "Omega_m = " << Omegam << endl;
  cerr << "Omega_L = " << OmegaL << endl;
  cerr << "bias factor = " << bg << endl;
  cerr << "Hubble constant = " << 100*h << "[km/s/Mpc]" << endl;
  cerr << "l_max = " << lmax << endl;
  /*
  cerr << "window function z0     = " << z0 << endl;
  cerr << "window function beta   = " << beta << endl;
  cerr << "window function lambda = " << lbda << endl;
  */
  
  ofstream output(outfname.c_str());

  // Initialize the spline over the matter 3D power spectrum  
  InitSpline(pkfname);

  // in case you've selected to work with the 2MASS/2MPZ photo-z distribuition as selection function,
  // a spline should also be performed for the photo-z histograms
  if (band>=6 && band<=10)
    InitSelectionFunctionSpline(band);

  // Start by writing in the output file, the monopole and dipole terms (zero by default)
  output << 0 << " " << 0. << endl; // monopole
  output << 1 << " " << 0. << endl; // dipole
  
  // Loop over the remaining multipoles (l<=10)
  
//File where time info of the loops will be output
//ofstream out2("outputs/z0_"+to_string(z0)+"_beta_"+to_string(beta)+"_lbda_"+to_string(lbda)+".txt");

//#pragma omp parallel for
  for (int l=48; l<=lmax; l++){
  
    InitWtSpline(l, OmegaL, Omegam, h);

    // Finally, get the cross-correlation
    double cttval = ctt(OmegaL, Omegam, l, h, bg);
    
    #pragma omp critical
    {
    cerr << l << " " << cttval << endl; // write cross-correlation to screen
    output << l << " " << cttval << endl; // write cross-correlation to output file
    }

    /*
    //Outputs time info to separate file
    out2 << "Loop for l=" << l << " started at " << start << endl;
    out2 << "Loop for l=" << l << " ended at " << stop << endl;
    */
  }

  // destroy splines
  DestroyWtSpline();
  DestroySpline();

  // close output file
  output.close();
  
  return 0;

}
