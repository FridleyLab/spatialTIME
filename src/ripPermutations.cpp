#include <Rcpp.h>
using namespace Rcpp;

// these are utility functions for calculating ripK permutations to get CSR distribution

// [[Rcpp::export]]
NumericMatrix cpp_matrix_subsetting(NumericMatrix m, NumericVector rows, NumericVector cols){
  
  int rl = rows.length();
  int cl = cols.length();
  NumericMatrix out(rl, cl);
  
  for (int i=0; i<cl; i++){
    NumericMatrix::Column org_c = m(_, cols[i]-1);
    NumericMatrix::Column new_c = out(_, i);
    for (int j=0; j<rl; j++){
      new_c[j] = org_c[rows[j]-1];
    }
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector compute_perms(NumericMatrix perms, double r_val, NumericMatrix distances, NumericMatrix edge, double area){
  int n = perms.cols();
  int cells = perms.rows();
  NumericVector K(n);
  
  for(int i=0; i<n; i++){
    NumericVector random_cells = perms.column(i);
    
    NumericMatrix subset_distances = cpp_matrix_subsetting(distances, random_cells, random_cells);
    NumericMatrix subset_edge = cpp_matrix_subsetting(edge, random_cells, random_cells);
    
    int rl = random_cells.length();
    int cl = rl;
    
    double summed_vals = 0.0;
    
    for(int j=0; j<cl; j++){
      for(int k=0; k<rl; k++){
        double val = subset_distances(j, k);
        if(val < r_val && val > 0){
          summed_vals =  summed_vals + subset_edge(j, k);
        }
      }
    }
    
    K[i] = (summed_vals * area)/(cells * (cells-1)); //(summed_vals * area)/(n * (n-1));
  }
  
  return(K);
}

// [[Rcpp::export]]
NumericVector compute_perms2(NumericMatrix perms, NumericVector r_range, NumericMatrix distances, NumericMatrix edge, double area){
  int n = perms.cols();
  int cells = perms.rows();
  NumericMatrix K(n, r_range.length());
  
  for(int r_index=0; r_index<r_range.length(); r_index++){
    double r_val = r_range[r_index];
    for(int i=0; i<n; i++){
      NumericVector random_cells = perms.column(i);
      
      NumericMatrix subset_distances = cpp_matrix_subsetting(distances, random_cells, random_cells);
      NumericMatrix subset_edge = cpp_matrix_subsetting(edge, random_cells, random_cells);
      
      int rl = random_cells.length();
      int cl = rl;
      
      double summed_vals = 0.0;
      
      for(int j=0; j<cl; j++){
        for(int k=0; k<rl; k++){
          double val = subset_distances(j, k);
          if(val < r_val && val > 0){
            summed_vals =  summed_vals + subset_edge(j, k);
          }
        }
      }
      
      K(i, r_index) = (summed_vals * area)/(cells * (cells-1));
    }
  }
  colnames(K) = r_range;
  return(K);
}



