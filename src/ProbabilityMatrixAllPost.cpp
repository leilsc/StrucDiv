// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
# include <Rcpp.h>
using namespace Rcpp ;



class NewProgressBar: public ProgressBar{
  public: // ====== LIFECYCLE =====
    
    /**
    * Main constructor
    */
    NewProgressBar()  { reset(); }
    
    ~NewProgressBar() {}
    
    public: // ===== main methods =====
      
      void display() {
        REprintf("Calculating gray level co-occurrence matrix\n");
        REprintf("0%%   10   20   30   40   50   60   70   80   90   100%%\n");
        REprintf("[----|----|----|----|----|----|----|----|----|----|\n");
        flush_console();
      }
      
      // will finalize display if needed
      void update(float progress) {
        _update_ticks_display(progress);
        if (_ticks_displayed >= _max_ticks)
          _finalize_display();
      }
      
      void end_display() {
        update(1);
        reset();
      }
      
      void reset() {
        _max_ticks = 50;
        _ticks_displayed = 0;
        _finalized = false;
      }
      
      
      protected: // ==== other instance methods =====
        
        // update the ticks display corresponding to progress
        void _update_ticks_display(float progress) {
          int nb_ticks = _compute_nb_ticks(progress);
          int delta = nb_ticks - _ticks_displayed;
          if (delta > 0) {
            _display_ticks(delta);
            _ticks_displayed = nb_ticks;
          }
          
        }
        
        void _finalize_display() {
          if (_finalized) return;
          
          REprintf("|\n");
          flush_console();
          _finalized = true;
        }
        
        int _compute_nb_ticks(float progress) {
          return int(progress * _max_ticks);
        }
        
        void _display_ticks(int nb) {
          for (int i = 0; i < nb; ++i) {
            REprintf("*");
            R_FlushConsole();
          }
        }
        
        // N.B: does nothing on windows
        void flush_console() {
          //if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
          R_FlushConsole();
          //endif
        }
        
private:
  int _max_ticks;   		// the total number of ticks to print
  int _ticks_displayed; 	// the nb of ticks already displayed
  bool _finalized;
  
};

// endif



// [[Rcpp::export(name = ".ProbabilityMatrixAllPost")]]
List ProbabilityMatrixAllPost(NumericMatrix vMat, NumericMatrix x, int d, 
                              int nrp, int nrp_big,
                              bool display_progress=true){
  
  List outList(vMat.nrow());
  
  NumericVector values_big = na_omit(x);
  NumericVector Values_big = sort_unique(values_big);
  
  NumericMatrix out_big(Values_big.length(), Values_big.length());
  
  for(int a = 0; a < x.nrow(); a++){
    for(int b = 0; b < x.ncol(); b++){
      for(int i = 0; i < out_big.nrow(); i++){
        for(int j = 0; j < out_big.ncol(); j++){
          
          if(b < x.ncol() - d){
            if(Values_big(i) == x(a,b) & Values_big(j) == x(a, b+d)){
              out_big(i,j) += 1;
            }
          }
          
          if( (b >= d) & (a < x.nrow() - d) ){
            if(Values_big(i) == x(a,b) & Values_big(j) == x(a+d, b-d)){
              out_big(i,j) += 1;
            }
          }
          
          if(a < x.nrow() - d){
            if(Values_big(i) == x(a,b) & Values_big(j) == x(a+d, b)){
              out_big(i,j) += 1;
            }
          }
          
          if( (b < x.ncol() - d) & (a < x.nrow() - d) ){
            if(Values_big(i) == x(a,b) & Values_big(j) == x(a+d, b+d)){
              out_big(i,j) += 1;
            }
          }
        }
      }
    }
  }
  
  NumericMatrix out2_big = transpose(out_big);
  
  for( int i = 0; i < out2_big.nrow(); i++ ){
    for( int j = 0; j < out2_big.ncol(); j++ ){
      
      // p.increment(); // update progress
      
      out2_big(i,j) += out_big(i,j);
      
    }
  }
  
  
  // NewProgressBar pb;
  // 
  // Progress p(vMat.nrow()*vMat.nrow(), display_progress, pb);
  

  for(int t = 0; t < vMat.nrow(); t++) {
    
    NumericVector rowvec = vMat(t,_);
    NumericMatrix xMat2( sqrt(vMat.ncol()), sqrt(vMat.ncol()), rowvec.begin() );
    NumericMatrix xMat = transpose(xMat2);
    NumericVector values = na_omit(rowvec);
    NumericVector Values = sort_unique(values);
    
    if(Values.length() == 0) continue;
    
    NumericMatrix out(Values.length(), Values.length());
    
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
        
        // p.increment(); // update progress
        
        out2(i,j) += out(i,j);

      }
    }
    
    colnames(out2) = Values;
    rownames(out2) = Values;
    
    
    IntegerVector ValPos = match(Values, Values_big);
    
    int rl = ValPos.length();
    NumericMatrix out2_big2(rl, rl);
    
    for (int i=0; i<rl; i++){
      NumericMatrix::Column org_c = out2_big(_, ValPos[i]-1);
      NumericMatrix::Column new_c = out2_big2(_, i);
      for (int j=0; j<rl; j++){
        new_c[j] = org_c[ValPos[j]-1];
      }
    }
    
    for( int i = 0; i < out2_big2.nrow(); i++ ){
      for( int j = 0; j < out2_big2.ncol(); j++ ){
        
        out2_big2(i,j) += out2(i,j);
        out2_big2(i,j) = out2_big2(i,j)/(nrp + nrp_big);
        
      }
    }
    
    
    colnames(out2_big2) = Values;
    rownames(out2_big2) = Values;
    
    outList[t] = out2_big2;
    
  }
  
  return outList;
  
}


