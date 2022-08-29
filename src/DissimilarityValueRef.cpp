# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export(name = ".DissimilarityValueRef")]]
NumericMatrix DissimilarityValueRef( NumericMatrix PMat, NumericVector xVal ){
  
  NumericMatrix out;
  
  CharacterVector xrows = rownames(PMat);
  CharacterVector xVals = rownames(PMat);
  std::transform(xrows.begin(), xrows.end(), xVals.begin(), std::atoi);
  
  NumericMatrix DisMat(PMat.nrow(), PMat.ncol());
  
  for (int m = 0; m < PMat.nrow(); m++) {
    for( int n = 0; n < PMat.ncol(); n++) {
      
      DisMat(m,n) = PMat(m,n) * fabs(xVal[m] - xVal[n]);
      
    }
  }
  
  out = DisMat;
  
  return(out);
  
}
