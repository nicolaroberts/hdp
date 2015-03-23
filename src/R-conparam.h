#ifndef _CONPARAM
#define _CONPARAM

#include "R-utils.h"

typedef struct {
  double alpha, alphaa, alphab;
  int numdp, *totalnd, *totalnt;
} CONPARAM;

CONPARAM *rReadConparamVector(SEXP cpin);

void rWriteConparamVector(SEXP result, CONPARAM *cpnew);

#endif