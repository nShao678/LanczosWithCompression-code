# LanczosWithCompression-code

This repo contains the implementation of Lanczos with compression proposed in [1] and scripts to reproduce the numerical experiments.

The main function of Lanczos with compression (Algorithm 2 in [1]) is LC.m.

Rational Krylov Toolbox (RKToolbox) is required, and can be downloaded at [http://guettel.com/rktoolbox/](http://guettel.com/rktoolbox/).

DFT matrices, including 'H2O','GaAsH6','SiO2','Si5H12','Ga10As10H30','Si34H36','Si41Ge41H72','Si87H76','Ge99H100', need to be downloaded from SuiteSparse Matrix Collection [https://sparse.tamu.edu/PARSEC](https://sparse.tamu.edu/PARSEC).
 

To reproduce numerical results, run main.m.


[1]: 
