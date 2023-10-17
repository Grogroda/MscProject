#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "chealpix.h"


int main(int argc, char **argv){

  if (argc<3){
    printf("Usage: ./apply_mask.exe [input map (fits)] [mask (fits)] [output map (fits)]\n");
    exit(1);
  }


  long int nside1, nside2; // nside for dust map
  char coordsys[10];
  char ordering[10];

  float *InputMap; // survey output map
  float *OutputMap; // survey output contrast map
  float *Mask;
 
  // loading input map
  InputMap = read_healpix_map(argv[1], &nside1, coordsys, ordering); 
  printf("Input map has nside = %ld \n",nside1);
  printf("Input map ordering: %s \n",ordering);

  // loading mask
  Mask = read_healpix_map(argv[2], &nside2, coordsys, ordering); 
  printf("Mask has nside = %ld \n",nside2);
  printf("Mask ordering: %s \n",ordering);


  long npix = nside2npix(nside1);

  OutputMap  = (float *)malloc(npix*sizeof(float));
  
  double healpix_bad_value=-1.6375e30;
  
  long ipnest;
  int i;
  double theta, phi;
  for (i=0; i<npix; i++){
    //pix2ang_nest(nside1, i, &theta, &phi);
    //ang2pix_nest(nside2, theta, phi, &ipnest);
    //ipnest = i;
    if (Mask[i]==1.)
      OutputMap[i] = InputMap[i];
    else
      OutputMap[i] = healpix_bad_value;
  }
  

  printf("Writing output map into healpix map format file %s \n",argv[3]);
  if (coordsys[0] == 'G'){
    printf("Galactic coordinates selected \n");
    write_healpix_map(OutputMap, nside1, argv[3], 0, "G");
  }
  else{
    printf("Equatorial coordinates selected \n");
    write_healpix_map(OutputMap, nside1, argv[3], 0, "C");
  }

  return 0;

}
