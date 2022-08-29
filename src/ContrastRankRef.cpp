# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export(name = ".ContrastRankRef")]]
NumericMatrix ContrastRankRef( NumericMatrix PMat ){
  
  NumericMatrix out;
  
  CharacterVector xrows = rownames(PMat);
  CharacterVector xVals = rownames(PMat);
  IntegerVector ValPos = match(xVals, xrows);
  std::transform(xrows.begin(), xrows.end(), xVals.begin(), std::atoi);
  
  NumericMatrix ConMat(PMat.nrow(), PMat.ncol());
  
  for (int m = 0; m < PMat.nrow(); m++) {
    for( int n = 0; n < PMat.ncol(); n++) {
      ConMat(m,n) = PMat(m,n) * ( pow(ValPos[m] - ValPos[n], 2) );
      
    }
  }
  
  out = ConMat;
  
  return(out);
  
}
