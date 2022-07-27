# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export(name = ".ProbabilityMatrixAll")]]
NumericMatrix ProbabilityMatrixAll(NumericMatrix xMat, int d, NumericVector Values){
  
  // allocate the matrix we will return
  NumericMatrix out(Values.length(), Values.length());
  
  int AllPairs = 2*( xMat.nrow()*(xMat.ncol() - d) + xMat.ncol()*(xMat.nrow() - d) +
                     2*(xMat.nrow() - d)*(xMat.ncol() - d) );
  
  for(int i = 0; i < out.nrow(); i++){
    for(int j = 0; j < out.ncol(); j++){
      for(int a = 0; a < xMat.nrow(); a++){
        for(int b = 0; b < xMat.ncol(); b++){
          
          if(b < xMat.ncol() - d){
            if(Values(i) == xMat(a,b) & Values(j) == xMat(a, b+d)){
              out(i,j) += 1;
            }
          }
          
          if( (b >= d) & (a < xMat.nrow() - d) ){
            if(Values(i) == xMat(a,b) & Values(j) == xMat(a+d, b-d)){
              out(i,j) += 1;
            }
          }
          
          if(a < xMat.nrow() - d){
            if(Values(i) == xMat(a,b) & Values(j) == xMat(a+d, b)){
              out(i,j) += 1;
            }
          }
          
          if( (b < xMat.ncol() - d) & (a < xMat.nrow() - d) ){
            if(Values(i) == xMat(a,b) & Values(j) == xMat(a+d, b+d)){
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
