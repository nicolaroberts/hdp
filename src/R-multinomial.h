#ifndef MULTINOMIAL
#define MULTINOMIAL

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>
#include "R-utils.h"


typedef struct {
  int numdim;
  double *eta;
} HH;

typedef int SS;

typedef int* QQ;


int numitems(QQ qq);

double marglikelihood(HH hh, QQ qq, SS ss);

double marglikelihoods(double *clik, HH hh, int numqq, QQ *qq, SS ss);

void adddata(HH hh, QQ qq, SS ss);

double adddatalik(HH hh, QQ qq, SS ss);

void deldata(HH hh, QQ qq, SS ss);

QQ newclass(HH hh);

SS *rReadSSVector(const SEXP ssdata);

SEXP rWriteSSVector(int numss, SS *ss);

HH rReadHH(const SEXP vector);

void FreeHH(HH hh);

QQ *rReadQQVector(HH hh, const SEXP qqmatrix, int maxnum);

SEXP rWriteQQVector(HH hh, int nn, int maxnum, QQ *qq);


#endif
