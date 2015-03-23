/*
 * Structure representing each Dirichlet process.  In particular, the
 * concentration parameter, the parent distribution, the beta weights, 
 * the numbers of data items/tables, and the statistics of data items 
 * directly associated with this DP.
 *
 * alpha        The concentration parameter itself.
 * numdata      The number of data items drawn directly from this DP.  This
 *              excludes virtual data items drawn as a result of data items
 *              associated with children DPs.
 * datacc       The cluster to which data items are associated with.  In
 *              matlab this is a row vector and datacc(ii) is the cluster of
 *              data item ii.
 * datass       The sufficient statistics of data items.  In matlab
 *              datass(:,ii) is the statistics for data item ii.
 * classnd      The number of data items in each component, including virtual
 *              data items.  There are numclass+1 entries, one for each
 *              component and an extra one for the empty class.  In matlab 
 *              this is a row vector, with cc'th entry being the number of 
 *              data items in component cc.
 * classnt      Same as classnd, except for number of tables.
 * beta         The weight associated with each component in this DP.  There
 *              are numclass+1 entries, one for each component, plus an
 *              additional one for the empty class.
 *
 * DP *mxReadDPVector(mxArray *mcell,int *dpstate,int mm)
 *              Reads in a DP struct from a matlab struct.  Allocates enough
 *              memory for mm components.
 * void mxWriteDPVector(mxArray *result,int nn,DP *dparray,int *dpstate,
 *                      int numclass)
 *              Writes a DP struct to a matlab struct.  Overwrites fields in
 *              result if necessary.  Frees memory allocated.
 */


#include "R-dp.h"

DP *rReadDPList(SEXP dpin, int *dpstate, int mm) {
  SEXP dpelt;
  DP *result, *dp;
  int ii, nn;
  nn = length(dpin);
  result = malloc(sizeof(DP)*nn);
  for ( ii = 0 ; ii < nn ; ii++ ) {
    dpelt       = PROTECT(VECTOR_ELT(dpin,ii));
    dp            = result+ii;
    dp->numdata   = asReal(rReadListElement(dpelt,"numdata"));
    dp->datass    = rReadSSVector(rReadListElement(dpelt,"datass"));
    if ( dpstate[ii] == ACTIVE || dpstate[ii] == FROZEN ) {
      //int cc;
      dp->classnd = rReadIntVector(rReadListElement(dpelt,"classnd"),mm,0,0);
      dp->classnt = rReadIntVector(rReadListElement(dpelt,"classnt"),mm,0,0);
      dp->datacc  = rReadIntVector(rReadListElement(dpelt,"datacc"),0,-1,0);
    } else {
      dp->classnd = NULL;
      dp->classnt = NULL;
      dp->datacc  = NULL;
    }
    if ( dpstate[ii] == ACTIVE ) {
      dp->alpha   = asReal(rReadListElement(dpelt,"alpha"));
      dp->beta    = rReadDoubleVector(rReadListElement(dpelt,"beta"),mm,0.0,0.0);
    } else {
      dp->alpha   = 0.0;
      dp->beta    = NULL;
    }
    UNPROTECT(1);
  }
  return result;
}

void rWriteDPList(SEXP result, int nn, DP *dpnew, int *dpstate, int numclass) {
  SEXP dpelt;
  DP *dp;
  int ii;
  rdebug2(3,"Write DP %d-vector, numclass=%d\n",nn,numclass);
  for ( ii = 0 ; ii < nn ; ii++ ) {
    dp = dpnew+ii;
    dpelt = PROTECT(VECTOR_ELT(result,ii));
    if ( dpstate[ii] == ACTIVE ) {
      rWriteListElement(dpelt,"alpha",rWriteRealScalar(dp->alpha));
      rWriteListElement(dpelt,"beta",rWriteDoubleVector(dp->beta, numclass, 0.0));
    } 

    if ( dpstate[ii] == ACTIVE || dpstate[ii] == FROZEN ) {
      //int cc;
      rWriteListElement(dpelt,"classnd", rWriteIntVector(dp->classnd,numclass,0));
      rWriteListElement(dpelt,"classnt", rWriteIntVector(dp->classnt,numclass,0));
      rWriteListElement(dpelt,"datacc", rWriteIntVector(dp->datacc,dp->numdata,1));
    }
  UNPROTECT(1);
  }  
  Free(dpnew);
}
