#include <gperftools/profiler.h>
#include <stdint.h>
#include <math.h>

#include "mcmc_utils.h"


namespace UtilFunctions {
  float fastlog2 (float x) {
    union { float f; uint32_t i; } vx = { x };
    union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
    float y = vx.i;
    y *= 1.1920928955078125e-7f;

    return y - 124.22551499f
            - 1.498030302f * mx.f 
            - 1.72587999f / (0.3520887068f + mx.f);
  }

  float fastlog (float x) {
    return 0.69314718f * fastlog2 (x);
  }

  int r_to_bool(SEXP x) {
    return Rcpp::as<bool>(x);
  };

  int r_to_int(SEXP x) {
    return Rcpp::as<int>(x);
  };

  double r_to_double(SEXP x) {
    return Rcpp::as<double>(x);
  };

  std::string r_to_string(SEXP x) {
    return Rcpp::as<std::string>(x);
  };

  std::vector<bool> r_to_vector_bool(SEXP x) {
    return Rcpp::as<std::vector<bool> >(x);
  };

  std::vector<int> r_to_vector_int(SEXP x) {
    return Rcpp::as<std::vector<int> >(x);
  };

  std::vector<double> r_to_vector_double(SEXP x) {
    return Rcpp::as<std::vector<double> >(x);
  };

  std::vector<std::string> r_to_vector_string(SEXP x) {
    return Rcpp::as<std::vector<std::string> >(x);
  };

  template <class T>
  std::vector<std::vector<T> > r_to_mat(Rcpp::List x) {
    int nrow = x.size();
    std::vector<std::vector<T> > x_mat(nrow);

    for(size_t i = 0; i < nrow; i++)
    {
      x_mat[i] = Rcpp::as<std::vector<T> > (x[i]);
    }
    
    return x_mat;
  };

  std::vector<std::vector<bool> > r_to_mat_bool(Rcpp::List x) {
    return r_to_mat<bool>(x);
  };

  std::vector<std::vector<int> > r_to_mat_int(Rcpp::List x) {
    return r_to_mat<int>(x);
  };

  std::vector<std::vector<double> > r_to_mat_double(Rcpp::List x) {
    return r_to_mat<double>(x);
  };

  template <class T>
  std::vector<std::vector<std::vector<T> > > r_to_array(Rcpp::List x) {
    int n_elements = x.size();
    UtilFunctions::print("Parsing Elements:", n_elements);
    std::vector<std::vector<std::vector<T> > > x_mat(n_elements);
    
    for(size_t i = 0; i < n_elements; i++) {
      Rcpp::List x_i(x[i]);
      int nrows = x_i.size();
      x_mat[i] = std::vector<std::vector<T> >(nrows);
      for(size_t j = 0; j < nrows; j++) {
        Rcpp::NumericVector x_i_j(x_i[j]);
        x_mat[i][j] = Rcpp::as<std::vector<T> >(x_i_j);
      }
    }
    
    return x_mat;
  };

  std::vector<std::vector<std::vector<bool> > > r_to_array_bool(Rcpp::List x) {
    return r_to_array<bool>(x);
  };

  std::vector<std::vector<std::vector<int> > > r_to_array_int(Rcpp::List x) {
    return r_to_array<int>(x);
  };

  std::vector<std::vector<std::vector<double> > > r_to_array_double(Rcpp::List x) {
    return r_to_array<double>(x);
  };

}

RcppExport SEXP start_profiler(SEXP str) {
  ProfilerStart(Rcpp::as<const char*>(str));
  return R_NilValue;
}

RcppExport SEXP stop_profiler() {
  ProfilerStop();
  return R_NilValue;
}