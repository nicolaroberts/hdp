#ifndef RANDUTILS
#define RANDUTILS

#include <R.h>
#include <math.h>
#include <stdlib.h>
#include "R-utils.h"

double randgamma(double rr);

int randnumtable(double alpha, int numdata);

void randdir(double *pi, double *alpha, int veclength, int skip);

double randbeta(double aa, double bb);

int randmult(double *pi, int veclength, int skip);

int randuniform(int numvalue);

double randconparam(double alpha, int numgroup, int *numdata, int *numtable,
                    double alphaa, double alphab, int numiter);




#endif
