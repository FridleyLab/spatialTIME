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

//' Compute the permutations of permuations of a distance matrix
//' 
//' @param perm1 a numeric matrix containing rows in `distances1` and `edge1` to use as positive cell locations
//' @param r_val a single r value which to count within for Ripley's K
//' @param distances1,edge1 matrices containing equal columns and rows for distances between cells and edge correction corresponding to cell pairs
//' @param area a single value that contains the area of the window around points
//' 
//' @return a vector of K statistics of length nrow(perms1)
//' @export
//' 
//' @examples
//' perms1 = as.matrix(sample(1:10, 5, replace = F))
//' r_val = 1
//' distances1 = matrix(nrow = 10, data = abs(rnorm(100)))
//' edge1 = distances1 #no edge correction
//' area = 100
//' K = compute_perms(perms1, r_val, distances1, edge1, area)
// [[Rcpp::export]] ////////// to make visible to package, put (rng = false) after explort before ]
NumericVector compute_perms(SEXP perms1, double r_val, SEXP distances1, SEXP edge1, double area){
  NumericMatrix perms(perms1), distances(distances1), edge(edge1), subset_distances, subset_edge;
  int n = perms.cols();
  int cells = perms.rows();
  NumericVector K(n);
  
  for(int i=0; i<n; i++){
    NumericVector random_cells = perms.column(i);
    
    subset_distances = cpp_matrix_subsetting(distances, random_cells, random_cells); //distance matrix for permutation
    subset_edge = cpp_matrix_subsetting(edge, random_cells, random_cells); //edge correction for permutation
    
    int rl = random_cells.length();
    int cl = rl;
    
    double summed_vals = 0.0;
    
    for(int j=0; j<cl; j++){ //column length
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
NumericVector compute_perms3(SEXP perms1, double r_val, SEXP distances1, SEXP edge1, double area){
  NumericMatrix perms(perms1), distances(distances1), edge(edge1), subset_distances, subset_edge;
  int n = perms.cols();
  int cells = perms.rows();
  NumericVector K(n), edge_c, val;
  
  for(int i=0; i<n; i++){
    NumericVector random_cells = perms.column(i);
    
    subset_distances = cpp_matrix_subsetting(distances, random_cells, random_cells); //distance matrix for permutation
    subset_edge = cpp_matrix_subsetting(edge, random_cells, random_cells); //edge correction for permutation
    
    int cl = random_cells.length();
    
    double summed_vals = 0.0;
    
    for(int j=0; j<cl; j++){ //column length
      NumericMatrix::Column dist_c = subset_distances(_, j);
      edge_c = subset_edge(_, j);
      val = edge_c[dist_c < r_val & dist_c > 0];
      summed_vals = summed_vals + sum(val);
      // for(int k=0; k<rl; k++){
      //   double val = subset_distances(j, k);
      //   if(val < r_val && val > 0){
      //     summed_vals =  summed_vals + subset_edge(j, k);
      //   }
      // }
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



