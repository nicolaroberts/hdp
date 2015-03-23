#ifndef _BASE
#define _BASE

#include <R.h>
#include <Rinternals.h>
#include "R-utils.h"
#include "R-multinomial.h"
#include <math.h>

typedef struct {
  int numclass, maxclass;
  HH hh;
  QQ *classqq;
  double *beta;
} BASE;

BASE *rReadBase(SEXP basein);

void rWriteBase(SEXP result, BASE *base);



#endif