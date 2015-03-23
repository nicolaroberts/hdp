#ifndef _HDP
#define _HDP

#include <math.h>
#include "R-utils.h"
#include "randutils.h"
#include "R-dp.h"
#include "R-base.h"
#include "R-conparam.h"

typedef struct {
  int numdp, numconparam;
  BASE *base;
  DP *dp; /* in top down order! */
  CONPARAM *conparam;
  double *clik;
  int *dpstate, *ppindex, *cpindex, *ttindex;
} HDP;


HDP *rReadHDP(SEXP hdpin);

void rWriteHDP(SEXP result, HDP *hdp);

int hdp_addclass(HDP *hdp);

int hdp_delclass(HDP *hdp, int cc);

double hdp_likelihood(HDP *hdp);

void hdp_randconparam(HDP *hdp, int numiter);

void hdp_randbeta(HDP *hdp, int jj);

void hdp_randclassnt(HDP *hdp, int jj);

void hdp_collecttotal(HDP *hdp, int jj);

void hdp_randdatacc(HDP *hdp, int jj);

void hdp_iterate(HDP *hdp, double *iterlik, int numiter, int doconparam, int dolik);

void hdp_dpactivate(HDP *hdp, int jj);

void hdp_dpholdout(HDP *hdp, int jj);

void hdp_predict(HDP *hdp, double *lik, int numburnin, int numsample, int numpredict, int *predictjj, int doconparam);

#endif

