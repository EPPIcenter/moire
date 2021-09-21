// #include <gperftools/profiler.h>
#include "mcmc_utils.h"

#include <R.h>
#include <Rcpp.h>
#include <vector>

namespace UtilFunctions
{
int r_to_bool(SEXP x) { return Rcpp::as<bool>(x); };

int r_to_int(SEXP x) { return Rcpp::as<int>(x); };

long int r_to_long_int(SEXP x) {return Rcpp::as<long int>(x); };

double r_to_double(SEXP x) { return Rcpp::as<double>(x); };

std::string r_to_string(SEXP x) { return Rcpp::as<std::string>(x); };

std::vector<bool> r_to_vector_bool(SEXP x)
{
    return Rcpp::as<std::vector<bool>>(x);
};

std::vector<int> r_to_vector_int(SEXP x)
{
    return Rcpp::as<std::vector<int>>(x);
};

std::vector<double> r_to_vector_double(SEXP x)
{
    return Rcpp::as<std::vector<double>>(x);
};

std::vector<std::string> r_to_vector_string(SEXP x)
{
    return Rcpp::as<std::vector<std::string>>(x);
};

std::vector<std::vector<bool>> r_to_mat_bool(
    Rcpp::Matrix<Rcpp::traits::r_sexptype_traits<bool>::rtype> x)
{
    return r_to_mat<bool>(x);
};

std::vector<std::vector<int>> r_to_mat_int(
    Rcpp::Matrix<Rcpp::traits::r_sexptype_traits<int>::rtype> x)
{
    return r_to_mat<int>(x);
};

std::vector<std::vector<double>> r_to_mat_double(
    Rcpp::Matrix<Rcpp::traits::r_sexptype_traits<double>::rtype> x)
{
    return r_to_mat<double>(x);
};

std::vector<std::vector<std::vector<bool>>> r_to_array_bool(Rcpp::List x)
{
    return r_to_array<bool>(x);
};

std::vector<std::vector<std::vector<int>>> r_to_array_int(Rcpp::List x)
{
    return r_to_array<int>(x);
};

std::vector<std::vector<std::vector<double>>> r_to_array_double(Rcpp::List x)
{
    return r_to_array<double>(x);
};

}  // namespace UtilFunctions

// RcppExport SEXP start_profiler(SEXP str) {
//   ProfilerStart(Rcpp::as<const char*>(str));
//   return R_NilValue;
// }
//
// RcppExport SEXP stop_profiler() {
//   ProfilerStop();
//   return R_NilValue;
// }
