# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export(name = ".ProbabilityMatrixVertical")]]
NumericMatrix ProbabilityMatrixVertical(NumericMatrix xMat, int d, NumericVector Values){
  
  int AllPairs = xMat.ncol() * 2 * (xMat.nrow() - d);
  
  // allocate the matrix we will return
  NumericMatrix out(Values.length(), Values.length());
  
  for(int i = 0; i < out.nrow(); i++){
    for(int j = 0; j < out.ncol(); j++){
      for(int a = 0; a < xMat.nrow(); a++){
        for(int b = 0; b < xMat.ncol(); b++){
          if(a < xMat.nrow() - d){
            
            if(Values(i) == xMat(a,b) & Values(j) == xMat(a+d, b)){
              out(i,j) += 1;
              
            }
          }
        }
      }
    }
  }
  
  NumericMatrix out2 = transpose(out);
  
  for( int i = 0; i < out2.nrow(); i++ ){
    for( int j = 0; j < out2.ncol(); j++ ){
      
      out2(i,j) += out(i,j);
      out2(i,j) = out2(i,j)/AllPairs;
      
    }
  }
  return out2;
  
  colnames(out2) = Values;
  rownames(out2) = Values;
  
}
