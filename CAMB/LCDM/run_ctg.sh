#!/bin/bash

#Change to the pk2ctg in the repo when quadrature integration is fixed

~/Dropbox/Estudos/Mestrado/Projeto/Cross-Correlation/code_xcross/para_arthur/bin/pk2ctg -pk matterPS.dat -wl 0.6847 -wm 0.3153 -h 0.6736 -lmax 96 -band 1 -bg 1.37 -out ctg2MASSband1_camb.dat 

~/Dropbox/Estudos/Mestrado/Projeto/Cross-Correlation/code_xcross/para_arthur/bin/pk2cgg -pk matterPS.dat -wl 0.6847 -wm 0.3153 -h 0.6736 -bg 1.37 -lmax 96 -band 1 -out cgg2MASSband1_camb.dat
