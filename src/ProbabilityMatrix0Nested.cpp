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



// [[Rcpp::export(name = ".ProbabilityMatrixHorizontalNested")]]
List ProbabilityMatrixHorizontalNested(NumericMatrix vMat, NumericMatrix vMat_big,
                                       int d, bool narm,
                                       bool display_progress=true){

  List outList(vMat.nrow());
  
  NewProgressBar pb;
  
  Progress p(vMat.nrow()*vMat.nrow(), display_progress, pb);
  
  int bb = sqrt(vMat.ncol());
  int bbb = sqrt(vMat_big.ncol());
  
  int AllPairs = ( bb*2*(bb - d) ) + ( bbb*2*(bbb - d) );
  

  for(int t = 0; t < vMat.nrow(); t++) {

    NumericVector rowvec = vMat(t,_);
    NumericMatrix xMat2( sqrt(vMat.ncol()), sqrt(vMat.ncol()), rowvec.begin() );
    NumericMatrix xMat = transpose(xMat2);  // transpose is fine because we
    // always have square matrices and we simply need to switch rows and columns

    NumericVector values = na_omit(rowvec);
    NumericVector Values = sort_unique(values);
    
    if(Values.length() == 0) continue;
    
    LogicalVector v = is_na(rowvec);
    
    if(narm==0 && any(v).is_true()) continue;

    // allocate the matrix we will return
    NumericMatrix out(Values.length(), Values.length());
    
    for(int a = 0; a < xMat.nrow(); a++){
      for(int b = 0; b < xMat.ncol(); b++){
        for(int i = 0; i < out.nrow(); i++){
          for(int j = 0; j < out.ncol(); j++){
            
            if(b < xMat.ncol() - d){

              if(Values(i) == xMat(a,b) & Values(j) == xMat(a, b+d)){
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
        
        p.increment(); // update progress
        
        out2(i,j) += out(i,j);
        out2(i,j) = out2(i,j);

      }
    }
    
    colnames(out2) = Values;
    rownames(out2) = Values;
  
  NumericVector rowvec_big = vMat_big(t,_);
  NumericMatrix xMat2_big( sqrt(vMat_big.ncol()), sqrt(vMat_big.ncol()), rowvec_big.begin() );
  NumericMatrix xMat_big = transpose(xMat2_big);  // transpose is fine because we
  // always have square matrices and we simply need to switch rows and columns
  
  NumericVector values_big = na_omit(rowvec_big);
  NumericVector Values_big = sort_unique(values_big);
  
  if(Values_big.length() == 0) continue;
  
  LogicalVector v_big = is_na(rowvec_big);
  
  if(narm==0 && any(v_big).is_true()) continue;
  
  // allocate the matrix we will return
  NumericMatrix out_big(Values_big.length(), Values_big.length());
  
  for(int a = 0; a < xMat.nrow(); a++){
    for(int b = 0; b < xMat.ncol(); b++){
      for(int i = 0; i < out.nrow(); i++){
        for(int j = 0; j < out.ncol(); j++){
          
          if(b < xMat.ncol() - d){
            
            if(Values(i) == xMat(a,b) & Values(j) == xMat(a, b+d)){
              out(i,j) += 1;
              
            }
          }
        }
      }
    }
  }
  
  NumericMatrix out2_big = transpose(out_big);
  
  for( int i = 0; i < out2_big.nrow(); i++ ){
    for( int j = 0; j < out2_big.ncol(); j++ ){
      
      out2_big(i,j) += out_big(i,j);
      
    }
  }
  
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
      
      p.increment(); // update progress
      
      out2_big2(i,j) += out2(i,j);
      out2_big2(i,j) = out2_big2(i,j)/AllPairs;
      
    }
  }
  
  colnames(out2_big2) = Values;
  rownames(out2_big2) = Values;
  
  outList[t] = out2_big2;
  
  }

  return outList;

}


