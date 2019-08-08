# thermal_iwf

This is the public repository for __thermal_iwf__. __thermal_iwf__ is the software accompanying the paper "Contemporary Meshfree Methods for Three Dimensional Heat Conduction Problems" published by Springer Nature in Archives of Computational Methods in Engineering. __thermal_iwf__ contains a wide array of meshless algorithms for approxmation of the laplacian in general and the solution of heat conduction problems in particlar. Namely:

* SPH
* SPH w/ Brookshaw Approximation
* CSPM
* RKPM
* PSE
* A new scheme presented by Fatehi et al. in the work "Fatehi, R., & Manzari, M. T. (2011). Error estimation in smoothed particle hydrodynamics and a new scheme for second derivatives. Computers and Mathematics with Applications, 61(2), 482–498"

__thermal_iwf__ is free software (GPLv3) and was written by Afrasiabi Mohamadreza (afrasiabi@ethz.ch) and Matthias Röthlin (mroethli@ethz.ch). __thermal_iwf__ is written in GNU C99 and does not require any dependencies. Makefiles for a debug and release build are provided. __thermal_iwf__ has been tested under Ubuntu 16.04, but does probably run on any flavor of Linux with minimal to no modification. 
