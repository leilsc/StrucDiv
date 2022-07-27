# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export(name = ".HomogeneityRankRef")]]
double HomogeneityRankRef( NumericMatrix PMat ){
  
  double out;
  
  CharacterVector xrows = rownames(PMat);
  CharacterVector xVals = rownames(PMat);
  IntegerVector ValPos = match(xVals, xrows);
  std::transform(xrows.begin(), xrows.end(), xVals.begin(), std::atoi);
  
  NumericMatrix HoMat(PMat.nrow(), PMat.ncol());
  
  for (int m = 0; m < PMat.nrow(); m++) {
    for( int n = 0; n < PMat.ncol(); n++) {
      HoMat(m,n) = PMat(m,n) / (1 + pow(ValPos[m] - ValPos[n], 2));
      
    }
  }
  
  out = sum(HoMat);
  
  return(out);
  
}
