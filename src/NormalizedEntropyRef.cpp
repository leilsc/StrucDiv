# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export(name = ".NormalizedEntropyRef")]]
NumericMatrix NormalizedEntropyRef( NumericMatrix PMat, NumericVector xVal, double nrp ){
  
  NumericMatrix out;
  
  double hmax = log(nrp);
  
  CharacterVector xrows = rownames(PMat);
  CharacterVector xVals = rownames(PMat);
  std::transform(xrows.begin(), xrows.end(), xVals.begin(), std::atoi);
  IntegerVector ValPos = match(xVals, xrows);
  
  NumericMatrix xPMatSub(xVal.length(), xVal.length());
  
  for (int s = 0; s < xPMatSub.ncol(); s++) {
    int pos = ValPos[s]-1;
    for( int t = 0; t < xPMatSub.nrow(); t++) {
      xPMatSub(t,s) = PMat(ValPos[t]-1, pos);
    }
  }
  
  NumericMatrix EntMat(PMat.nrow(), PMat.ncol());
  
  for (int m = 0; m < xPMatSub.nrow(); m++) {
    for( int n = 0; n < xPMatSub.ncol(); n++) {
      
      if(xPMatSub(m,n) == 0) EntMat(m,n) = xPMatSub(m,n) * 0;
      else EntMat(m,n) = xPMatSub(m,n) * (-log(xPMatSub(m,n))) / hmax;
      
    }
  }
  
  out = EntMat;
  
  return(out);
  
}
