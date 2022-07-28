# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export(name = ".HomogeneityValueRef")]]
double HomogeneityValueRef( NumericMatrix PMat, NumericVector xVal ){
  
  double out;
  
  CharacterVector xrows = rownames(PMat);
  CharacterVector xVals = rownames(PMat);
  std::transform(xrows.begin(), xrows.end(), xVals.begin(), std::atoi);
  
  NumericMatrix HoMat(PMat.nrow(), PMat.ncol());
  
  for (int m = 0; m < PMat.nrow(); m++) {
    for( int n = 0; n < PMat.ncol(); n++) {
      HoMat(m,n) = PMat(m,n) / (1 + pow(xVal[m] - xVal[n], 2));
      
    }
  }
  
  out = sum(HoMat);
  
  return(out);
  
}
