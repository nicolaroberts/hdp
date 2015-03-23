#ifndef RUTILS
#define RUTILS

#include <math.h>
#include <R.h>
#include <Rinternals.h>

SEXP rReadListElement(const SEXP list, const char *str);
void rWriteListElement(SEXP list, const char *str, SEXP newelement);
SEXP rWriteRealScalar(double var);
SEXP rWriteIntScalar(int var);
int *rReadIntVector(SEXP rvec, int number, int shift, int init);
double *rReadDoubleVector(SEXP rvec, int number, double shift, double init);
SEXP rWriteIntVector(int *var, int len, int shift);
SEXP rWriteDoubleVector(double *var, int len, double shift);

#define max(x1,x2) ( (x1) < (x2) ? (x2) : (x1) )

#ifndef NODEBUG
int DEBUG;
#define rdebug0(num,string)       if(DEBUG>=(num))Rprintf(string);
#define rdebug1(num,string,a)     if(DEBUG>=(num))Rprintf(string,a);
#define rdebug2(num,string,a,b)   if(DEBUG>=(num))Rprintf(string,a,b);
#define rdebug3(num,string,a,b,c) if(DEBUG>=(num))Rprintf(string,a,b,c);
#define rdebug4(num,string,a,b,c,d) if(DEBUG>=(num))Rprintf(string,a,b,c,d);
#define rdebugarray(num,string,str,array,length) { \
  if (DEBUG >= (num)) { \
    int ii; \
    Rprintf("%s: ",string); \
    for ( ii = 0 ; ii < length ; ii++) { \
      Rprintf(str,array[ii]); \
      if (ii % 5 == 4) { \
        Rprintf(" : "); \
      } else { \
        Rprintf(" ");  \
      } \
    } \
    Rprintf("\n"); \
  } \
}
#else
#define rdebug0(num,string)     
#define rdebug1(num,string,a)     
#define rdebug2(num,string,a,b)  
#define rdebug3(num,string,a,b,c)
#define rdebug4(num,string,a,b,c,d)
#define rdebugarray(num,string,str,array,length) 
#endif

#endif