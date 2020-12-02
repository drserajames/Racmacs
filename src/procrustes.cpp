
#include <RcppArmadillo.h>
#include "procrustes.h"
using namespace Rcpp;

// Define a procrustes transformation
// [[Rcpp::export]]
Procrustes ac_procrustes(
  arma::mat X,
  arma::mat Xstar,
  bool translation,
  bool dilation
){

  if(X.n_rows != Xstar.n_rows){ Rcpp::stop("X and Xstar do not have same number of rows."); }
  if(X.n_cols != Xstar.n_cols){ Rcpp::stop("X and Xstar do not have same number of columns."); }

  int n = X.n_rows;
  int m = X.n_cols;

  arma::mat J = arma::mat(n,n);
  J.eye();
  if(translation){
    J -= (1.0/n);
  }

  arma::mat tX = arma::trans(X);
  arma::mat tXstar = arma::trans(Xstar);

  arma::mat C = tXstar * J * X;

  arma::vec svd_d;
  arma::mat svd_u;
  arma::mat svd_v;
  arma::svd(
    svd_u,
    svd_d,
    svd_v,
    C
  );

  arma::mat R = svd_v * svd_u;
  double s = 1.0;

  if(dilation){

    arma::mat mat1 = tXstar * J * X * R;
    arma::mat mat2 = tX * J * X;

    double s_numer = 0.0;
    double s_denom = 0.0;

    for(int i=0; i<m; i++){
      s_numer += mat1(i,i);
      s_denom += mat2(i,i);
    }

    s = s_numer / s_denom;

  }

  arma::mat tt = arma::mat(m, 1, arma::fill::zeros);
  if(translation){

    arma::mat mat1 = arma::mat(n, 1, arma::fill::ones);
    arma::mat tmatmid = arma::trans(Xstar - s * X * R);
    tt = ((1.0/n) * tmatmid)*mat1;

  }

  Procrustes out;
  out.R = R;
  out.tt = tt;
  out.s = s;
  return out;

}

// Apply a procrustes transformation
arma::mat ac_apply_procrustes(
    arma::mat coords,
    Procrustes p
){

  return transform_coords(
    coords,
    p.R,
    p.tt,
    p.s
  );

}

// Align coordinates via procrustes
// [[Rcpp::export]]
arma::mat ac_align_coords(
    arma::mat source,
    arma::mat target,
    bool translation = true,
    bool dilation = false
){

  Procrustes p = ac_procrustes(
    source,
    target,
    translation,
    dilation
  );

  return ac_apply_procrustes(source, p);

}

// Apply a coordinate transformation
arma::mat transform_coords(
  const arma::mat coords,
  const arma::mat rotation,
  const arma::mat translation,
  const double scaling
) {

  // Work out maximum dims
  int dims = arma::max(
    arma::uvec{
      coords.n_cols,
      rotation.n_cols,
      translation.n_rows
    }
  );

  // Expand matrices to match maximum dimensions
  arma::mat tcoords(coords.n_rows, dims);
  tcoords.cols(0, coords.n_cols - 1) = coords;

  arma::mat trotation(dims, dims, arma::fill::eye);
  if(rotation.n_rows > 0){
    trotation.submat(
      0, 0,
      rotation.n_rows - 1, rotation.n_cols - 1
    ) = rotation;
  }

  arma::mat ttranslation(coords.n_rows, dims, arma::fill::zeros);
  for(int i=0; i<ttranslation.n_rows; i++){
    for(int j=0; j<translation.n_rows; j++){
      ttranslation(i,j) = translation(j,0);
    }
  }

  // Perform the transformation
  arma::mat out = (scaling*tcoords)*trotation + ttranslation;
  return out;

}

