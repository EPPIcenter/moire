#pragma once

#ifndef MCMC_UTILS_H_
#define MCMC_UTILS_H_

#include <Rcpp.h>

#define OVERFLO 1e100
#define UNDERFLO 1e-100

namespace UtilFunctions
{
float fastlog2(float x);

float fastlog(float x);

int r_to_bool(SEXP x);

int r_to_int(SEXP x);

double r_to_double(SEXP x);

std::string r_to_string(SEXP x);

std::vector<bool> r_to_vector_bool(SEXP x);

std::vector<int> r_to_vector_int(SEXP x);

std::vector<double> r_to_vector_double(SEXP x);

std::vector<std::string> r_to_vector_string(SEXP x);

template <class T>
std::vector<std::vector<T>> r_to_mat(Rcpp::List x);

std::vector<std::vector<bool>> r_to_mat_bool(Rcpp::List x);

std::vector<std::vector<int>> r_to_mat_int(Rcpp::List x);

std::vector<std::vector<double>> r_to_mat_double(Rcpp::List x);

template <class T>
std::vector<std::vector<std::vector<T>>> r_to_array(Rcpp::List x);

std::vector<std::vector<std::vector<bool>>> r_to_array_bool(Rcpp::List x);

std::vector<std::vector<std::vector<int>>> r_to_array_int(Rcpp::List x);

std::vector<std::vector<std::vector<double>>> r_to_array_double(Rcpp::List x);

template <class T>
void print(T x)
{
    Rcpp::Rcout << x << "\n";
    R_FlushConsole();
};

template <class T, class... Args>
void print(T first, Args... args)
{
    Rcpp::Rcout << first << " ";
    print(args...);
}

template <class T>
void print_vector(std::vector<T> v)
{
    Rcpp::Rcout << '[';
    for (size_t i = 0; i < v.size() - 1; i++)
    {
        Rcpp::Rcout << v[i] << ", ";
    }
    Rcpp::Rcout << v[v.size() - 1] << "]"
                << "\n";
}
}  // namespace UtilFunctions

#endif  // MCMC_UTILS_H_
