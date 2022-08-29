# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export(name = ".DissimilarityRankRef")]]
NumericMatrix DissimilarityRankRef( NumericMatrix PMat ){
  
  NumericMatrix out;
  
  CharacterVector xrows = rownames(PMat);
  CharacterVector xVals = rownames(PMat);
  IntegerVector ValPos = match(xVals, xrows);
  std::transform(xrows.begin(), xrows.end(), xVals.begin(), std::atoi);
  
  NumericMatrix DisMat(PMat.nrow(), PMat.ncol());
  
  for (int m = 0; m < PMat.nrow(); m++) {
    for( int n = 0; n < PMat.ncol(); n++) {
      DisMat(m,n) = PMat(m,n) * abs(ValPos[m] - ValPos[n]);
      
    }
  }
  
  out = DisMat;
  
  return(out);
  
}
