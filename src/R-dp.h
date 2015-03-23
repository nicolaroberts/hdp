#ifndef _DP
#define _DP

#include "R-utils.h"
#include "R-multinomial.h"

typedef struct DP{
  double alpha, *beta;
  int *classnd, *classnt;
  int numdata, *datacc;
  SS *datass;
} DP;


DP *rReadDPList(SEXP dpin, int *dpstate, int mm);
void rWriteDPList(SEXP result, int nn, DP *dpnew, int *dpstate, int numclass);


#define ACTIVE 2
#define FROZEN 1
#define HELDOUT 0

#endif