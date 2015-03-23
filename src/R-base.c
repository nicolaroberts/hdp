/*
 * Structure representing the base distribution, in particular, the 
 * parameters for the base distribution, and samples from the base 
 * distribution (one sample per class).
 *
 * numclass     The number of samples from the base distribution (number of
 *              clusters in HDP mixture).
 * maxclass     The number of samples for which we have allocated memory.
 *              This is initialized to (numclass+2)*2.
 * hh           The parameters for the base distribution.
 * classqq      The sufficient statistics for data items associated with
 *              each cluster.  In matlab, classqq(:,kk) refers to the
 *              statistics for class kk.  There is numclass+1 represented
 *              clusters (one additional cluster with no data associated 
 *              with it).
 * beta         Only used internally.  A vector of numclass 0's following by
 *              one 1.  In matlab this is a row vector.
 *
 * BASE *rReadBase(rArray *basein);
 *              Reads in a BASE struct from a matlab struct.
 * void rWriteBase(rArray *result,BASE *base);
 *              Writes a BASE struct to a matlab struct.  Overwrites ListElements
 *              in result if necessary.  Frees memory allocated.
 */


#include "R-base.h"


BASE *rReadBase(SEXP basein) {
  BASE *result;
  int ii, maxclass;
  result = malloc(sizeof(BASE));
  result->numclass = asInteger(rReadListElement(basein,"numclass")); 
  result->maxclass = maxclass = (result->numclass+2) * 2;
  result->hh       = rReadHH(rReadListElement(basein,"hh"));
  result->classqq  = rReadQQVector(result->hh,rReadListElement(basein,"classqq"), result->maxclass);
  result->beta = malloc(sizeof(double)*maxclass);
  for ( ii = 0 ; ii < maxclass ; ii++) 
    result->beta[ii] = 0.0;
  result->beta[result->numclass] = 1.0;
  return result;
}


void rWriteBase(SEXP result, BASE *base) {
  rWriteListElement(result,"numclass",rWriteIntScalar(base->numclass));
  rWriteListElement(result,"classqq",rWriteQQVector(base->hh,base->numclass+1,base->maxclass,base->classqq));
  FreeHH(base->hh);
  Free(base->beta);
  Free(base);
}


