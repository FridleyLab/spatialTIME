#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector getWindow(int x, int size) {
  
  if(x < size){
    NumericMatrix m(1,2);
    m(0,0) = 1;
    m(0,1) = x;
    return(m);
  }
  
  double quo = x / size;
  double rem = x % size;
  
  NumericVector start (quo+1, 0);
  start[0] = 1;
  NumericVector end (quo+1, 0);
  
  for(int i = 0; i < quo; i++){
    end[i] = (i+1)*size;
  }
  
  for(int i = 1; i <= quo; i++){
    start[i] = i*size+1;
  }
  
  if(rem == 0){
    start.erase(quo);
    end.erase(quo);
  }
  if(rem > 0){
    end[quo] = rem + (quo * size);
  }
  
  NumericMatrix m(start.length(), 2);
  
  m(_,0) = start;
  m(_,1) = end;
  
  return(m);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
slide(42, 100)
slide(200, 100)
*/
