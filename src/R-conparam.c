/*
 * Structure representing a concentration parameter in the HDP.  In
 * particular, we represent the parameter itself, and parameters for
 * a gamma prior over it, the number of DPs using the concentration
 * parameter, and the total number of data items and tables for each DP.
 *
 * alpha        The concentration parameter itself.  Has to be nonnegative.
 * alphaa       The shape parameter for the gamma prior over alpha.
 * alphab       The inverse scale parameter for the gamma prior over alpha.
 * numdp        The number of DPs using this concentration parameter.
 * totalnd      The number of data items in each DP.  In matlab this is a row
 *              vector with jj'th entry belonging to the jj'th DP associated
 *              this concentration parameter.
 * totalnt      Same as totalnd for number of tables in each DP.
 *
 * CONPARAM *mxReadConparamVector(mxArray *mcell);
 *              reads in a CONPARAM struct from a matlab struct. 
 * void mxWriteConparamVector(mxArray *result,CONPARAM *cparray);
 *              writes a CONPARAM struct to a matlab struct.  Overwrites 
 *              fields if necessary.  Frees memory allocated.
 */


#include "R-conparam.h"

CONPARAM *rReadConparamVector(SEXP cpin) {
  SEXP cplistelt;
  CONPARAM *result, *cp;
  int ii, nn;
  nn = length(cpin);
  result = malloc(sizeof(CONPARAM)*nn);
  for ( ii = 0 ; ii < nn ; ii++ ) {
    cplistelt = PROTECT(VECTOR_ELT(cpin,ii));
    cp  = result + ii;
    cp->alpha = asReal(rReadListElement(cplistelt,"alpha"));
    cp->alphaa = asReal(rReadListElement(cplistelt,"alphaa"));
    cp->alphab = asReal(rReadListElement(cplistelt,"alphab"));
    cp->numdp  = asInteger(rReadListElement(cplistelt,"numdp"));
    cp->totalnd = rReadIntVector(rReadListElement(cplistelt,"totalnd"),0,0,0);
    cp->totalnt = rReadIntVector(rReadListElement(cplistelt,"totalnt"),0,0,0);
    UNPROTECT(1);
  }
  return result;
}

void rWriteConparamVector(SEXP result, CONPARAM *cpnew) {
  SEXP cplistelt;
  CONPARAM *cp;
  int ii, nn;
  nn = length(result);
  for ( ii = 0 ; ii < nn ; ii++ ) {
    cp = cpnew + ii;
    cplistelt = PROTECT(VECTOR_ELT(result,ii));
    rWriteListElement(cplistelt,"alpha",rWriteRealScalar(cp->alpha));
    rWriteListElement(cplistelt,"totalnd",rWriteIntVector(cp->totalnd,cp->numdp,0));
    rWriteListElement(cplistelt,"totalnt",rWriteIntVector(cp->totalnt,cp->numdp,0));
    UNPROTECT(1);
  }
  Free(cpnew);
}
