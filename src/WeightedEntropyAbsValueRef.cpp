# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export(name = ".WeightedEntropyAbsValueRef")]]
double WeightedEntropyAbsValueRef( NumericMatrix PMat ){
  
  double out;
  
  CharacterVector xrows = rownames(PMat);
  CharacterVector xVals = rownames(PMat);
  NumericVector xVal = rownames(PMat);
  std::transform(xrows.begin(), xrows.end(), xVals.begin(), std::atoi);
  
  NumericMatrix EntMat(PMat.nrow(), PMat.ncol());
  
  for (int m = 0; m < PMat.nrow(); m++) {
    for( int n = 0; n < PMat.ncol(); n++) {
      if(PMat(m,n) == 0) EntMat(m,n) = PMat(m,n) * 0;
      else EntMat(m,n) = std::abs( xVal[m] - xVal[n] )*( PMat(m,n) * (-log(PMat(m,n))) );
      
    }
  }
  
  out = sum(EntMat);
  
  return(out);
  
}
