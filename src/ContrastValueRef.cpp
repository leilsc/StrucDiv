# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export(name = ".ContrastValueRef")]]
double ContrastValueRef( NumericMatrix PMat ){
  
  double out;
  
  CharacterVector xrows = rownames(PMat);
  CharacterVector xVals = rownames(PMat);
  NumericVector xVal = rownames(PMat);
  std::transform(xrows.begin(), xrows.end(), xVals.begin(), std::atoi);
  
  NumericMatrix ConMat(PMat.nrow(), PMat.ncol());
  
  for (int m = 0; m < PMat.nrow(); m++) {
    for( int n = 0; n < PMat.ncol(); n++) {
      ConMat(m,n) = PMat(m,n) * ( pow(xVal[m] - xVal[n], 2) );
      
    }
  }
  
  out = sum(ConMat);
  
  return(out);
  
}
