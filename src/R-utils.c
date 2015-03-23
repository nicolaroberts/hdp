
#include "R-utils.h"

SEXP rReadListElement(const SEXP list, const char *str) {
  SEXP elmt = R_NilValue;
  SEXP names = getAttrib(list, R_NamesSymbol);
  int i; 
  for ( i = 0; i < length(list); i++ )
    if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) { 
      elmt = VECTOR_ELT(list, i); 
      break; 
    }

  if ( elmt == R_NilValue )
    error("%s missing from list", str);

  if (DEBUG>=3) Rprintf("Read %s.\n",str); 
  return elmt;
}

void rWriteListElement(SEXP list, const char *str, SEXP newelement) {
  SEXP names = getAttrib(list, R_NamesSymbol);
  int i; 
  for ( i = 0; i < length(list); i++ )
    if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) { 
      SET_VECTOR_ELT(list, i, newelement);
      if (DEBUG>=3) Rprintf("Write %s.\n",str);
    }
}

SEXP rWriteRealScalar(double var) { 
  SEXP result = PROTECT(ScalarReal(var));
  UNPROTECT(1);
  return result;
}

SEXP rWriteIntScalar(int var) { 
  SEXP result = PROTECT(ScalarInteger(var));
  UNPROTECT(1);
  return result;
}

int *rReadIntVector(SEXP rvec, int number, int shift, int init){
  int *copy;
  int *result;
  int ii;
  number = max(number,length(rvec));
  result = malloc(sizeof(int)*number);
  copy = INTEGER(rvec);
  for ( ii = 0 ; ii < length(rvec) ; ii++ )
      result[ii] = copy[ii] + shift;
  for ( ii = length(rvec) ; ii < number ; ii++ )
      result[ii] = init;
  if (DEBUG>=3) {
    Rprintf("Value = ");
    for ( ii = 0 ; ii < number ; ii++ )
      Rprintf("%d ",result[ii]);
    Rprintf("\n");
    }
  return result;
}

double *rReadDoubleVector(SEXP rvec, int number, double shift, double init){
  double *copy;
  double *result;
  int ii;
  number = max(number,length(rvec));
  result = malloc(sizeof(double)*number);
  copy = REAL(rvec);
  for ( ii = 0 ; ii < length(rvec) ; ii++ )
      result[ii] = copy[ii] + shift;
  for ( ii = length(rvec) ; ii < number ; ii++ )
      result[ii] = init;
  if (DEBUG>=3) {
    Rprintf("Value = ");
    for ( ii = 0 ; ii < number ; ii++ )
      Rprintf("%g ",result[ii]);
    Rprintf("\n");
    }
  return result;
}

SEXP rWriteIntVector(int *var, int len, int shift) { 
  SEXP result = PROTECT(allocVector(INTSXP, len));
  int *copy = INTEGER(result);
  int ii;
  for ( ii = 0 ; ii < len ; ii++ )
    copy[ii] = var[ii] + shift;
  if (DEBUG>=3) {
    Rprintf("Value = ");
    for ( ii = 0 ; ii < len ; ii++ )
      Rprintf("%d ",var[ii]);
    Rprintf("\n");
  }
  Free(var);
  UNPROTECT(1);
  return result;
}

SEXP rWriteDoubleVector(double *var, int len, double shift) { 
  SEXP result = PROTECT(allocVector(REALSXP, len));
  double *copy = REAL(result);
  int ii;
  for ( ii = 0 ; ii < len ; ii++ )
    copy[ii] = var[ii] + shift;
  if (DEBUG>=3) {
    Rprintf("Value = ");
    for ( ii = 0 ; ii < len ; ii++ )
      Rprintf("%g ",var[ii]);
    Rprintf("\n");
  }
  Free(var);
  UNPROTECT(1);
  return result;
}
