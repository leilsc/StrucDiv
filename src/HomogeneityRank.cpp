// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
# include <Rcpp.h>
using namespace Rcpp ;



class MyProgressBar: public ProgressBar{
  public: // ====== LIFECYCLE =====
    
    /**
     * Main constructor
     */
    MyProgressBar()  { reset(); }
    
    ~MyProgressBar() {}
    
    public: // ===== main methods =====
      
      void display() {
        REprintf("Calculating horizontal structural diversity\n");
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



// [[Rcpp::export(name = ".HomogeneityRank")]]
NumericVector HomogeneityRank( NumericMatrix Hetx, List PMat,
                               bool  narm,
                               bool display_progress=true){
  
  NumericVector out(Hetx.nrow());
  
  MyProgressBar pb;
  
  Progress p(Hetx.nrow()*Hetx.nrow(), display_progress, pb);
  
  for(int i = 0; i < Hetx.nrow(); i++){
    
    NumericVector x = Hetx(i,_);
    
    LogicalVector v = is_na(x);
    
    if(narm==0 && any(v).is_true()) {
      
      out[i] = NA_REAL;
      
    }
    
    else {
      NumericVector xVal_ = na_omit(x);
      if(xVal_.length()==0) {
        
        out[i] = NA_REAL;
        
      }
      
      else {
        NumericVector xVal = sort_unique(xVal_);
        NumericMatrix xPMat(xVal.length(), xVal.length());
        rownames(xPMat) = xVal;
        CharacterVector xVals = rownames(xPMat);
        NumericMatrix PMat_ = PMat[i];
        CharacterVector xrows = rownames(PMat_);
        IntegerVector ValPos = match(xVals, xrows);
        
        NumericMatrix xPMatSub(xVal.length(), xVal.length());
        
        for (int s = 0; s < xPMatSub.ncol(); s++) {
          int pos = ValPos[s]-1;
          for( int t = 0; t < xPMatSub.nrow(); t++) {
            xPMatSub(t,s) = PMat_(ValPos[t]-1, pos);
          }
        }
        
        NumericMatrix HoMat(ValPos.length(), ValPos.length());
        
        for (int m = 0; m < xPMatSub.nrow(); m++) {
          for( int n = 0; n < xPMatSub.ncol(); n++) {
            
            p.increment(); // update progress
            
            HoMat(m,n) = xPMatSub(m,n) / (1 + pow(ValPos[m] - ValPos[n], 2));
          }
        }
        
        out[i] = sum(HoMat);
        
      }
    }
  }
  
  return(out);
  
}
