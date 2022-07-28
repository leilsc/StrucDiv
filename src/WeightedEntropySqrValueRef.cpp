# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export(name = ".WeightedEntropySqrValueRef")]]
double WeightedEntropySqrValueRef( NumericMatrix PMat, NumericVector xVal ){
  
  double out;
  
  CharacterVector xrows = rownames(PMat);
  CharacterVector xVals = rownames(PMat);
  std::transform(xrows.begin(), xrows.end(), xVals.begin(), std::atoi);
  
  NumericMatrix EntMat(PMat.nrow(), PMat.ncol());
  
  for (int m = 0; m < PMat.nrow(); m++) {
    for( int n = 0; n < PMat.ncol(); n++) {
      if(PMat(m,n) == 0) EntMat(m,n) = PMat(m,n) * 0;
      else EntMat(m,n) = pow( xVal[m] - xVal[n], 2 )*( PMat(m,n) * (-log(PMat(m,n))) );
      
    }
  }
  
  out = sum(EntMat);
  
  return(out);
  
}
