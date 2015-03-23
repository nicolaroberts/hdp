
#include "R-multinomial.h"

int numitems(QQ qq) {
  return qq[0];
}

double marglikelihood(HH hh, QQ qq, SS ss) {
  return log( (hh.eta[ss]+qq[ss])/(hh.eta[0]+qq[0]) );
}

double marglikelihoods(double *clik, HH hh, int numqq, QQ *qq, SS ss) {
  double etas, eta0;
  int ii;
  eta0 = hh.eta[0];
  etas = hh.eta[ss];
  for ( ii = 0 ; ii < numqq ; ii++ ) 
    clik[ii] = (etas + qq[ii][ss]) / (eta0 + qq[ii][0]);
  return 0;
}

void adddata(HH hh, QQ qq, SS ss) {
  qq[ss] ++;
  qq[0] ++;
}

double adddatalik(HH hh, QQ qq, SS ss) {
  return log( (hh.eta[ss] + (qq[ss] ++)) / (hh.eta[0] + (qq[0] ++)) );
}

void deldata(HH hh, QQ qq, SS ss) {
  qq[ss] --;
  qq[0] --;
}

QQ newclass(HH hh) {
  int ii;
  QQ result;
  result = malloc(sizeof(int)*(hh.numdim+1));
  for ( ii = 0 ; ii <= hh.numdim ; ii++ ) 
    result[ii] = 0;
  return result;
}

SS *rReadSSVector(const SEXP ssdata) {
  SS *result;
  int *pr;
  int ii;
  pr = INTEGER(ssdata);
  result = malloc(sizeof(SS)*length(ssdata));
  for ( ii = 0 ; ii < length(ssdata) ; ii++ )
    result[ii] = pr[ii];
  return result;
}

SEXP rWriteSSVector(int numss, SS *ss) {
  SEXP result;
  int *pr;
  int ii;
  result = PROTECT(allocVector(INTSXP,numss));
  pr = INTEGER(result);
  for ( ii = 0 ; ii < numss ; ii++ )
    pr[ii] = ss[ii];
  Free(ss);
  return result;
}

HH rReadHH(const SEXP vector) {
  double *pr, sum;
  int ii;
  HH result;
  result.numdim = length(vector);
  result.eta    = malloc(sizeof(double)*(1+result.numdim));
  sum           = 0.0;
  pr            = REAL(vector);
  for ( ii = 0 ; ii < result.numdim ; ii++ ) 
    sum += result.eta[ii+1]  = pr[ii];
  result.eta[0] = sum;
  return result;
}

void FreeHH(HH hh) {
  Free(hh.eta);
}

QQ *rReadQQVector(HH hh, const SEXP qqmatrix, int maxnum) {
  int *pr;
  int ii, jj, mm, nn, sum;
  QQ *result;
  SEXP dims = PROTECT(getAttrib(qqmatrix, R_DimSymbol));
  mm = INTEGER(dims)[0];
  nn = INTEGER(dims)[1];
  if ( mm != hh.numdim ) error("Number of dimensions don't match.");
  maxnum = max(maxnum,nn);
  result = malloc(sizeof(QQ)*maxnum);

  pr = INTEGER(qqmatrix); 

  for ( jj = 0 ; jj < nn ; jj++ ) {
    result[jj] = newclass(hh);
    sum = 0;
    for ( ii = 0 ; ii < mm ; ii++ ) {
      sum += result[jj][ii+1] = pr[ii+jj*mm];
    }
    result[jj][0] = sum;
  }
  for ( jj = nn ; jj < maxnum ; jj++ ) 
    result[jj] = newclass(hh);
  UNPROTECT(1);
  return result;
}

SEXP rWriteQQVector(HH hh, int nn, int maxnum, QQ *qq) {
  SEXP result;
  int *pr;
  int ii,jj;
  result = PROTECT(allocMatrix(INTSXP,hh.numdim,nn));
  pr = INTEGER(result);
  rdebug2(3,"Write qq (%dx%d): ",hh.numdim,nn);
  for ( jj = 0 ; jj < nn ; jj ++ ) {
    for ( ii = 0 ; ii < hh.numdim ; ii++ )
      pr[ii+jj*hh.numdim] = qq[jj][ii+1];
    Free(qq[jj]);
    rdebug1(4,"%d ",jj);
  }
  rdebug0(3,"Free qq: ");
  for ( jj = nn ; jj < maxnum ; jj ++ ) {
    Free(qq[jj]);
    rdebug1(4,"%d ",jj);
  }
  Free(qq);
  UNPROTECT(1);
  rdebug0(3,"\n");
  return result;
}
