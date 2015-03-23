/*
 * Structure for hierarchical Dirichlet process mixtures.  Contains a base
 * distribution, concentration parameters, DPs arranged in a tree, and
 * data items for each DP.
 *
 * numdp        The number of DPs.
 * numconparam  The number of concentration parameters.
 * base         The base distribution.
 * conparam     The list of concentration parameters.
 * dp           The list of DPs.  This is arranged in topological order
 *              consistent with the tree structure (first DP is root).
 * dpstate      List of states of DPs.
 * ppindex      List of indices of parents of DPs.  If this is 0 in matlab or
 *              -1 in C, this means the base distribution is the parent.
 * cpindex      List of indices of concentration parameters used by DPs
 * ttindex      List of indices into totalnd/nt arrays used by conparams.
 * clik         Internal use, vector of doubles.
 *
 * hdp *rReadHDP(const rArray *mcell)
 *              Reads in a HDP from a matlab struct.
 * void rWriteHDP(rArray *result, HDP *hdp)
 *              Writes a HDP to a matlab struct.  Overwrites fields in
 *              result if necessary.  Frees memory allocated.
 * int hdp_addclass(HDP *hdp)
 *              Adds one new class, increment maxclass if necessary, 
 *              returns new number of classes.
 * int hdp_delclass(HDP *hdp,int cc)
 *              Deletes class cc, returns new number of classes. 
 * double hdp_likelihood(HDP *hdp)
 *              Calculates log probability of generating observations in 
 *              active DPs, given current hidden state configuration. 
 * void hdp_randconparam(HDP *hdp,int numiter)
 *              Samples concentration parameters. 
 * void hdp_randbeta(HDP *hdp,int jj)
 *              Sample beta (mixing proportions) variables.
 * void hdp_randclassnt(HDP *hdp,int jj)
 *              Sample classnt (number of tables) variables.
 * void hdp_collecttotal(HDP *hdp,int jj)
 *              Sum number of tables and data items per DP. 
 * void hdp_dpactivate(HDP *hdp, int jj)
 *              Activates DP jj.  Assumes parent is already activates.  
 *              Otherwise crashes. 
 * void hdp_dpholdout(HDP *hdp, int jj)
 *              Holds out DP jj.  Assumes parent is activated.  
 *              Otherwise crashes.
 * void hdp_iterate(HDP *hdp, double *iterlik, int numiter, int doconparam,
 *     int dolik)
 *              Samples HDP using beta auxiliary variable method for numiter
 *              iterations.  If doconparam>0 then also samples concentration
 *              parameters, for doconparam number of iterations each time.
 *              If dolik is 1 then also calculates and returns the log
 *              probability of data items associated with activated DPs in 
 *              iterlik.
 * void hdp_predict(HDP *hdp, double *lik, int numburnin, int numsample, 
 *     int numpredict, int *predictjj, int doconparam)
 *              Estimates the log probability of data items in DPs given by
 *              indices in predictjj.  There are numpredict such DPs.  For
 *              each such DP we individually estimate the log probability of 
 *              the data items in that DP using the harmonic average
 *              technique of Kass and Raftery 1995, using numburnin number 
 *              of burn in iterations, followed by numsample iterations where
 *              we collect the log probabilities.  Returns result in lik,
 *              which in matlab is a numsample x numpredict array of doubles.
 */

#include "R-hdp.h"


/***************************************************************************/
HDP *rReadHDP(SEXP hdpin) {
  HDP *result;
  result = malloc(sizeof(HDP));
  rdebug0(1,"Reading HDP.\n");
  result->numdp       = asInteger(rReadListElement(hdpin,"numdp"));
  result->numconparam = asInteger(rReadListElement(hdpin,"numconparam"));
  result->base        = rReadBase(rReadListElement(hdpin,"base"));
  result->dpstate     = rReadIntVector(rReadListElement(hdpin,"dpstate"),result->numdp,0,0);
  result->ppindex     = rReadIntVector(rReadListElement(hdpin,"ppindex"),result->numdp,-1,0);
  result->cpindex     = rReadIntVector(rReadListElement(hdpin,"cpindex"),result->numdp,-1,0);
  result->ttindex     = rReadIntVector(rReadListElement(hdpin,"ttindex"),result->numdp,-1,0);
  result->dp          = rReadDPList(rReadListElement(hdpin,"dp"),
                            result->dpstate,result->base->maxclass);
  result->conparam    = rReadConparamVector(rReadListElement(hdpin,"conparam"));
  result->clik        = malloc(sizeof(double)*result->base->maxclass);
  return result;
}

void rWriteHDP(SEXP result, HDP *hdp) {
  rdebug0(1,"Writing HDP.\n");
  rWriteDPList(rReadListElement(result,"dp"),hdp->numdp,hdp->dp,hdp->dpstate,hdp->base->numclass+1);
  rWriteConparamVector(rReadListElement(result,"conparam"),hdp->conparam);
  rWriteBase(rReadListElement(result,"base"),hdp->base);
  rWriteListElement(result,"dpstate",rWriteIntVector(hdp->dpstate,hdp->numdp,0));
  //Free(hdp->ppindex);
  //Free(hdp->cpindex);
  //Free(hdp->ttindex);
  //Free(hdp->clik);
  //Free(hdp);
}

/***************************************************************************/
#define hdp_extendclass(pointer,start,size,type,zero) { \
  type *pp, *pe; \
  pointer = realloc(pointer, sizeof(type)*size); \
  pp = pointer + start; \
  pe = pointer + size; \
  while ( pp < pe ) { *pp = zero; pp++; } \
} 
int hdp_addclass(HDP *hdp) {
  BASE *base;
  DP *alldp, *dp;
  HH hh;
  int numclass, maxclass, numdp, pp, *dpstate, *ppindex;
  double *beta, alpha;

  int jj;
  double bb,b1,b2;

  base     = hdp->base;
  alldp    = hdp->dp;
  numdp    = hdp->numdp;
  numclass = base->numclass;
  maxclass = base->maxclass;
  dpstate  = hdp->dpstate;
  ppindex  = hdp->ppindex;
  rdebug1(1,"addclass %d. ",numclass);

  /* Stick-break beta weights */
  for ( jj = 0 ; jj < numdp ; jj++ ) {
    dp     = alldp + jj;
    if ( dpstate[jj] == ACTIVE ) {
      pp     = ppindex[jj];
      alpha  = dp->alpha;
      if ( pp == -1 ) {
        b1   = randgamma(1.0);
        b2   = randgamma(alpha);
      } else {
        beta = alldp[pp].beta;
        b1   = randgamma(alpha*beta[numclass]);
        b2   = randgamma(alpha*beta[numclass+1]);
      }
      beta   = dp->beta;
      if ((b1+b2)==0){
        beta[numclass] = 0;
        beta[numclass+1] = 0;
      } else {
        bb     = beta[numclass] / (b1 + b2);
        beta[numclass]   = bb * b1;
        beta[numclass+1] = bb * b2;
      }
      rdebug3(3,"\n  jj %d: beta %1.3g %1.3g.",jj,
          beta[numclass],beta[numclass+1]);
    } else {
      rdebug2(3,"\n  jj %d: skipped. state %d",jj,dpstate[jj]);
    }
  }
  rdebug0(3,"\n");
  base->beta[numclass]   = 0.0;
  base->beta[numclass+1] = 1.0;
  base->numclass         = numclass += 1;

  if (numclass+1 >= maxclass) {
    /* increase pool if necessary */
    base->maxclass = maxclass *= 2;
    rdebug1(1,"maxclass %d. ",maxclass);
    numdp = hdp->numdp;
    hh    = base->hh;
    hdp_extendclass(base->classqq,numclass+1,maxclass,QQ,newclass(hh));
    hdp_extendclass(hdp->clik,    numclass+1,maxclass,double,0.0);
    hdp_extendclass(base->beta,   numclass+1,maxclass,double,0.0);
    for ( jj = 0 ; jj < numdp ; jj++ ) {
      dp    = alldp + jj;
      if ( dpstate[jj] == ACTIVE ) {
        hdp_extendclass(dp->beta,   numclass+1,maxclass,double,0.0);
      }
      if ( dpstate[jj] != HELDOUT ) {
        hdp_extendclass(dp->classnd,numclass+1,maxclass,int,0);
        hdp_extendclass(dp->classnt,numclass+1,maxclass,int,0);
      }
    }
  }
 return numclass;
} 


/***************************************************************************/
#define hdp_deleteclass(array,start,size,type) { \
  type *pp, *pe, tmpvar; \
  pp     = array + start; \
  pe     = array + size; \
  tmpvar = *pp; \
  while ( pp < pe ) { *pp = *(pp+1); pp++; } \
  *pe = tmpvar; \
}
int hdp_delclass(HDP *hdp, int cc) {
  BASE *base;
  DP *alldp, *dp;
  int numclass, numdata, numdp, *datacc, *dpstate; 

  int ii, jj;

  rdebug1(1,"delclass %d. ",cc);
  base      = hdp->base;
  alldp     = hdp->dp;
  numclass  = base->numclass;
  numdp     = hdp->numdp;
  dpstate   = hdp->dpstate;
  hdp_deleteclass(base->classqq,cc,numclass,QQ);
  for ( jj = 0 ; jj < numdp ; jj++ ) {
    dp      = alldp + jj;
    if ( dpstate[jj] == ACTIVE ) {
      hdp_deleteclass(dp->beta,   cc,numclass,double);
    }
    if ( dpstate[jj] != HELDOUT ) {
      rdebugarray(3,"\n  old classnd","%2d",dp->classnd,numclass+1);
      hdp_deleteclass(dp->classnd,cc,numclass,int);
      rdebugarray(3,"  new classnd","%2d",dp->classnd,numclass+1);
      rdebugarray(3,"  old classnt","%2d",dp->classnt,numclass+1);
      hdp_deleteclass(dp->classnt,cc,numclass,int);
      rdebugarray(3,"  new classnt","%2d",dp->classnt,numclass+1);
      numdata = dp->numdata;
      datacc  = dp->datacc;
      for ( ii = 0 ; ii < numdata ; ii++ ) 
        datacc[ii] -= ( datacc[ii] > cc );
    }
  }
  base->beta[numclass]       = 0.0;
  base->beta[numclass-1]     = 1.0;
  base->numclass = numclass -= 1;
  return numclass;
}

/***************************************************************************/
double hdp_likelihood(HDP *hdp) {
  DP *alldp, *dp;
  SS *datass;
  HH hh;
  QQ *classqq;
  int *datacc, numdp, numdata, *dpstate; 
  
  int ii, jj;
  double lik;

  rdebug0(1,"likelihood: ");
  /* remove all data items from their classes */
  numdp   = hdp->numdp;
  alldp   = hdp->dp;
  dpstate = hdp->dpstate;
  hh      = hdp->base->hh;
  classqq = hdp->base->classqq;
  for ( jj = 0 ; jj < numdp ; jj++ ) {
    dp = alldp + jj;
    if ( dpstate[jj] == ACTIVE ) {
      numdata = dp->numdata;
      datass  = dp->datass;
      datacc  = dp->datacc;
      for ( ii = 0 ; ii < numdata ; ii++ )
        deldata(hh,classqq[datacc[ii]],datass[ii]);
    }
  }
  /* calculate likelihood and add data items back in */
  lik = 0.0;
  for ( jj = 0 ; jj < numdp ; jj++ ) {
    dp = alldp + jj;
    if ( dpstate[jj] == ACTIVE ) {
      numdata = dp->numdata;
      datass  = dp->datass;
      datacc  = dp->datacc;
      for ( ii = 0 ; ii < numdata ; ii++ ) 
        lik += adddatalik(hh,classqq[datacc[ii]],datass[ii]);
    }
  }
  rdebug1(1,"%1.3g. ",lik);
  return lik;
} /* hdp_likelihood */

/***************************************************************************/
void hdp_randconparam(HDP *hdp, int numiter) {
  DP *alldp, *dp;
  CONPARAM *allconparam, *conparam;
  int numconparam, numdp, *dpstate, *cpindex; 

  int cp, jj;

  rdebug0(2,"sample conparam. "); 

  /* sample alpha */
  allconparam = hdp->conparam;
  numconparam = hdp->numconparam;
  dpstate     = hdp->dpstate;
  cpindex     = hdp->cpindex;

  rdebug0(3,"\n");
  for ( cp = 0 ; cp < numconparam ; cp++ ) {
    conparam = allconparam + cp;
    rdebug2(3,"  cp %d: oldalpha %1.3g ",cp,conparam->alpha);
    conparam->alpha = randconparam(conparam->alpha,
        conparam->numdp,conparam->totalnd,conparam->totalnt,
        conparam->alphaa,conparam->alphab,numiter);
    rdebug1(3," newalpha %1.3g.\n",conparam->alpha); 
  }

  /* update DP's alphas */
  alldp       = hdp->dp;
  numdp       = hdp->numdp;
  for ( jj = 0 ; jj < numdp ; jj++ ) {
    dp = alldp + jj;
    if ( dpstate[jj] == ACTIVE ) {
      cp        = cpindex[jj];
      dp->alpha = allconparam[cp].alpha;
    }
  }
} 

/***************************************************************************/
void hdp_randbeta(HDP *hdp, int jj) {
  BASE *base;
  DP *alldp, *dp;
  int *classnd, numclass;
  double alpha, *beta;

  int cc, pp;
  double *clik;

  base     = hdp->base;
  numclass = base->numclass;
  alldp    = hdp->dp;
  dp       = alldp + jj;
  clik     = hdp->clik;

  rdebug0(2,"sample beta. ");
  classnd = dp->classnd;
  alpha   = dp->alpha;
  pp      = hdp->ppindex[jj];
  beta = ( pp == -1 ) ? base->beta : alldp[pp].beta;
  for ( cc = 0 ; cc <= numclass ; cc++ )
    clik[cc] = classnd[cc] + alpha*beta[cc];
  randdir(dp->beta,clik,numclass+1,1);
  rdebugarray(3,"\n  beta","%1.3g",dp->beta,numclass+1); 
}

void hdp_randclassnt(HDP *hdp, int jj) {
  DP *alldp, *dp;
  CONPARAM *conparam;
  int numclass, numdp;
  double alpha, *beta;
  int *classnd, *classnt, *pclassnd;

  /* temp variables */ 
  int cc, pp; // deleted tt, nd, sp, nt because compiler said they were unused

  conparam = hdp->conparam;
  numclass = hdp->base->numclass;
  numdp    = hdp->numdp;
  alldp    = hdp->dp;
  dp       = alldp + jj;

  pp       = hdp->ppindex[jj];
  alpha    = dp->alpha;
  classnd  = dp->classnd;
  classnt  = dp->classnt;

  rdebug0(2,"sample classnt. ");
  if ( pp == -1 ) {
    for ( cc = 0 ; cc < numclass ; cc++ ) {
      classnt[cc] = ( classnd[cc] > 0 );
    }
    rdebugarray(3,"\n  new classnt","%d",classnt,numclass+1);
  } else {
    pclassnd = alldp[pp].classnd;
    beta     = alldp[pp].beta;
    for ( cc = 0 ; cc < numclass ; cc++ ) {
      pclassnd[cc] -= classnt[cc];
      pclassnd[cc] += ( classnt[cc] = 
                        randnumtable(alpha*beta[cc],classnd[cc]) );
    }
    rdebugarray(3,"  new classnt","%d",classnt,numclass+1);
  } 
} 

void hdp_collecttotal(HDP *hdp, int jj) {
  DP *dp;
  CONPARAM *conparam;
  int numclass, *classnd, *classnt, cp, tt;
  int cc, nd, nt;

  conparam = hdp->conparam;
  numclass = hdp->base->numclass;
  dp       = hdp->dp + jj;
  
  cp       = hdp->cpindex[jj];
  tt       = hdp->ttindex[jj];
  classnd  = dp->classnd;
  classnt  = dp->classnt;
  nd = 0;
  nt = 0;
  
  for ( cc = 0 ; cc < numclass ; cc++ ) {
    nd += classnd[cc];
    nt += classnt[cc];
  }
  conparam[cp].totalnd[tt] = nd;
  conparam[cp].totalnt[tt] = nt;
} 

void hdp_randdatacc(HDP *hdp, int jj) {
  /* variables to be read in from hdp struct */
  BASE *base;
  DP *alldp, *dp; //deleted unused variable *par
  HH hh;
  QQ *classqq;
  SS *datass, ss;
  int numclass, *datacc, numdp, numdata, pp;
  double alpha, *beta;
  int *classnd, *classnt;

  /* temp variables */ 
  int ii, cc; //deleted unused variable kk
  double *clik;

  base        = hdp->base;
  numclass    = base->numclass;
  hh          = base->hh;
  classqq     = base->classqq;
  numdp       = hdp->numdp;
  alldp       = hdp->dp;
  dp          = alldp + jj;
  clik        = hdp->clik;

  rdebug1(2,"jj %d\n",jj);
  rdebug0(2,"Sample datacc. ");
  pp         = hdp->ppindex[jj];
  numdata    = dp->numdata;
  datass     = dp->datass;
  datacc     = dp->datacc;
  classnd    = dp->classnd;
  classnt    = dp->classnt;
  alpha      = dp->alpha;
  beta       = ( pp == -1 ) ? base->beta : alldp[pp].beta;

  for ( ii = 0 ; ii < numdata ; ii++ ) {
    ss = datass[ii];
    cc = datacc[ii];
    rdebug4(3,"\n  DP %d item %d oldcc %d numclass %d.\n",
        jj,ii,cc,numclass);

    /* sample class assignment for each data item */ 
    rdebug0(3,"  remove data.\n");
    deldata(hh,classqq[cc],ss);
    rdebugarray(3,"  old classnd","%d",classnd,numclass+1);
    classnd[cc] -= 1;

    marglikelihoods(clik,hh,numclass+1,classqq,ss);
    for ( cc = 0 ; cc <= numclass ; cc++ ) 
      clik[cc] *= classnd[cc] + alpha*beta[cc];

    rdebugarray(3,"  clik","%1.3g",clik,numclass+1);
    datacc[ii] = cc = randmult(clik,numclass+1,1);

    rdebug1(3,"add data. newcc %d. ",cc);
    adddata(hh,classqq[cc],ss);
    classnd[cc] += 1;
    rdebugarray(3,"new classnd","%d",classnd,numclass+1);

    /* increase number of classes if necessary */
    if ( cc == numclass ) {
      numclass = hdp_addclass(hdp);
      classqq  = base->classqq;
      classnd  = dp->classnd;
      classnt  = dp->classnt;
      beta     = pp == -1 ? base->beta : alldp[pp].beta;
      clik     = hdp->clik;
    }
  } 
} 

/***************************************************************************/
void hdp_iterate(HDP *hdp, double *iterlik,
    int numiter, int doconparam, int dolik) {
  BASE *base;
  QQ *classqq;
  DP *alldp;
  int numdp, *dpstate;
  int jj, cc, iter;

  numdp   = hdp->numdp;
  alldp   = hdp->dp;
  dpstate = hdp->dpstate;

  for ( iter = 0 ; iter < numiter ; iter++ ) {
    rdebug1(1,"iter %d: ",iter);
    for ( jj = numdp-1 ; jj >= 0 ; jj-- ) {
      rdebug2(2,"\n DP %d state %d: ",jj,dpstate[jj]);
      if ( dpstate[jj] == ACTIVE ) {
        hdp_randdatacc(hdp, jj);
        hdp_randclassnt(hdp, jj);
        hdp_collecttotal(hdp, jj);
      }
    }
    rdebug0(2,"\n");
    for ( jj = 0 ; jj < numdp ; jj++ ) {
      if ( dpstate[jj] == ACTIVE ) hdp_randbeta(hdp, jj);
    }

    /* delete empty classes, only after randclassnt, randbeta for consistency */
    base = hdp->base;
    classqq = base->classqq;
    for ( cc = base->numclass-1 ; cc >= 0 ; cc-- ) {
      if ( numitems(classqq[cc]) == 0 ) hdp_delclass(hdp,cc);
    }

    if ( doconparam > 0 ) hdp_randconparam(hdp,doconparam);

    if ( dolik == 1 ) iterlik[iter] = hdp_likelihood(hdp);

    rdebug0(1,"\n");
  } /* iter */
}

void hdp_dpactivate(HDP *hdp, int jj) {
  BASE *base;
  CONPARAM *conparam;
  DP *alldp, *dp;
  HH hh;
  QQ *classqq;
  double *beta;
  int numclass, maxclass, pp, cp, *classnd, *classnt, numdata, *datacc, 
        *dpstate;
  SS *datass;
  int cc, ii;

  base     = hdp->base;
  maxclass = base->maxclass;
  alldp    = hdp->dp;
  dpstate  = hdp->dpstate;
  dp       = alldp + jj;

  rdebug1(1,"Activating DP %d.\n",jj);
  if ( dpstate[jj] == HELDOUT ) {
    dp->classnd = classnd = malloc(sizeof(int)*maxclass);
    dp->classnt = classnt = malloc(sizeof(int)*maxclass);
    for ( cc = 0 ; cc < maxclass ; cc++ ) {
      classnd[cc] = 0;
    }

    numdata  = dp->numdata;
    numclass = base->numclass;
    hh       = base->hh;
    classqq  = base->classqq;
    datass   = dp->datass;
    dp->datacc = datacc = malloc(sizeof(int)*numdata);
    for ( ii = 0 ; ii < numdata ; ii++ ) {
      datacc[ii] = cc = randuniform(numclass);
      rdebug1(1,"cc %d\n",cc);
      adddata(hh,classqq[cc],datass[ii]);
      classnd[cc] += 1;
    }
    for ( cc = 0 ; cc < maxclass ; cc++ ) 
      classnt[cc] = classnd[cc];

    pp = hdp->ppindex[jj];
    if ( pp > -1 ) {
      classnd = alldp[pp].classnd;
      for ( cc = 0 ; cc < numclass ; cc++ )
        classnd[cc] += classnt[cc];
      hdp_collecttotal(hdp,pp);
    }
  }

  conparam  = hdp->conparam;
  cp        = hdp->cpindex[jj];
  dp->alpha = conparam[cp].alpha;
  dp->beta  = beta = malloc(sizeof(double)*maxclass);
  for ( cc = 0 ; cc < maxclass ; cc++ ) 
    beta[cc] = 0.0;
  hdp_randbeta(hdp,jj);
  hdp_collecttotal(hdp,jj);
  dpstate[jj] = ACTIVE;
}

void hdp_dpholdout(HDP *hdp, int jj) {
  CONPARAM *conparam;
  BASE *base;
  DP *alldp, *dp;
  HH hh;
  QQ *classqq;
  SS *datass;
  int numclass, numdata, pp, cp, tt, *classnd, *classnt, *pclassnd, *datacc,
        *dpstate, *ppindex;

  int cc, ii; 

  rdebug1(1,"Holding out DP %d.\n",jj);
  alldp    = hdp->dp;
  ppindex  = hdp->ppindex;
  dpstate  = hdp->dpstate;

  base     = hdp->base;
  hh       = base->hh;
  classqq  = base->classqq;

  dp       = alldp + jj;
  numdata  = dp->numdata;
  datacc   = dp->datacc;
  datass   = dp->datass;
  for ( ii = 0 ; ii < numdata ; ii++ )
    deldata(hh,classqq[datacc[ii]],datass[ii]);

  numclass = base->numclass;
  pp       = ppindex[jj];
  classnd  = dp->classnd;
  classnt  = dp->classnt;
  if ( pp > -1 ) {
    pclassnd = alldp[pp].classnd;
    for ( cc = 0 ; cc < numclass ; cc++ ) {
      pclassnd[cc] -= classnt[cc];
    }
    while ( pp > -1 ) {
      hdp_randclassnt(hdp,pp);
      hdp_collecttotal(hdp,pp);
      pp = ppindex[pp];
    }
  }

  conparam = hdp->conparam;
  cp       = hdp->cpindex[jj];
  tt       = hdp->ttindex[jj];
  conparam[cp].totalnd[tt] = 0;
  conparam[cp].totalnt[tt] = 0;

  Free(dp->classnd);
  Free(dp->classnt);
  if ( dpstate[jj] == ACTIVE )
    Free(dp->beta);

  dpstate[jj] = HELDOUT;
}

void hdp_predict(HDP *hdp, double *lik, int numburnin, int numsample, 
    int numpredict, int *predictjj, int doconparam) {
  int jj;


  for ( jj = 0 ; jj < numpredict ; jj++ )
    hdp_dpholdout(hdp,predictjj[jj]);
  for ( jj = 0 ; jj < numpredict ; jj++ ) {
    hdp_dpactivate(hdp,predictjj[jj]);
    hdp_iterate(hdp, NULL, numburnin, doconparam, 0);
    hdp_iterate(hdp, lik, numsample, doconparam, 1);
    lik += numsample;
    hdp_dpholdout(hdp,predictjj[jj]);
  }
  for ( jj = 0 ; jj < numpredict ; jj++ )
    hdp_dpactivate(hdp,predictjj[jj]);

}

