#!/bin/bash

# Compiling the object files:

g++ -c survey.cc -o survey.o -I/usr/local/include -fPIC
g++ -c cosmology.cc -o cosmology.o -I/usr/local/include -fPIC
g++ -c correlations.cc -o correlations.o -I/usr/local/include -fPIC 

# Creating library:

g++ -shared -o libcorrelations.so correlations.o cosmology.o survey.o -L/usr/local/lib -lm -lgsl -lgslcblas
