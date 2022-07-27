# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export(name = ".DissimilarityValueRef")]]
double DissimilarityValueRef( NumericMatrix PMat ){
  
  double out;
  
  CharacterVector xrows = rownames(PMat);
  CharacterVector xVals = rownames(PMat);
  NumericVector xVal = rownames(PMat);
  std::transform(xrows.begin(), xrows.end(), xVals.begin(), std::atoi);
  
  NumericMatrix DisMat(PMat.nrow(), PMat.ncol());
  
  for (int m = 0; m < PMat.nrow(); m++) {
    for( int n = 0; n < PMat.ncol(); n++) {
      
      DisMat(m,n) = PMat(m,n) * fabs(xVal[m] - xVal[n]);
      
    }
  }
  
  out = sum(DisMat);
  
  return(out);
  
}
