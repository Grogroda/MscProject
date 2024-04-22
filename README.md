# Cosmological Constraints with the Sachs-Wolfe Integrated Effect
In this repository I'll attempt to organize and make my Master's project available to my collaborators (in an organized way). It's still a work in progress, there is much to be organized before it can be effectively used, but once it does it would be great to only work in these codes through this repository. To work with this repo, one should read my Master's project, which includes many references necessary to understand what many of the codes present here do.

## General info

I'll not track test files, test scripts, some data files and I'll also not track the codes I write when I'm learning a new module, library or code in general. I'll also not track bin directories. All directory paths will be relative paths because I'm assuming everyone who gets access to these codes will maintain this repo's structure, which makes everything simpler. In case you need to transfer some of this to a cluster for instance, or if you just downloaded a few of the codes, it's likely you will need to adapt the code to account for this. 

For Linux users, I recommend using the "tree" software to take a look at the directory tree of this repo before working on it. Compilation and running instructions are available in every main directory of the repo, all of them taking into account I use Ubuntu 22.04 as my OS with GCC 11.3.0.  

## About the main directories

- cross-correlation: Contains codes that calculate the CMB-galaxy map cross correlation and does many other jobs using these results, such as minimizing the likelihood for a set of parameters of the selection function;
- CAMB: Contains codes meant to use the [CAMB](https://camb.info/) module to generate power spectra according to desired models and following certain directives, many of the files created with these codes may be used in the calculation of cross-correlations.

## Software Needed

- [CAMB](https://camb.info/): Recommended normal installation using pip: pip3 install camb;
- [Cobaya](https://cobaya.readthedocs.io/en/latest/installation.html): Recommended normal installation with MPI;
- [HEALPix](https://healpix.sourceforge.io/): More specifically, I've used its Python version called [healpy](https://healpy.readthedocs.io/en/latest/);
- [GSL](https://www.gnu.org/software/gsl/): Recommended standard installation (let the script choose the installation path, probably using /usr/local/lib and /usr/local/include);
- Common Python modules/libraries: This project uses [Numpy](https://numpy.org/), [Matplotlib](https://matplotlib.org/) and [Pandas](https://pandas.pydata.org/) for various numerical and graphical tasks.
