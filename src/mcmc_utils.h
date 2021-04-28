#pragma once

#ifndef MCMC_UTILS_H_
#define MCMC_UTILS_H_

#include <Rcpp.h>

#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
#include <Rinterface.h>
#endif

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
std::vector<std::vector<T>> r_to_mat(
    Rcpp::Matrix<Rcpp::traits::r_sexptype_traits<T>::rtype> x)
{
    int nrow = x.nrow();
    int ncol = x.ncol();
    std::vector<std::vector<T>> x_mat(nrow);

    for (size_t i = 0; i < nrow; i++)
    {
        for (size_t j = 0; j < ncol; j++)
        {
            x_mat[i].push_back(x.at(i, j));
        }
    }

    return x_mat;
};

std::vector<std::vector<bool>> r_to_mat_bool(
    Rcpp::Matrix<Rcpp::traits::r_sexptype_traits<bool>::rtype> x);

std::vector<std::vector<int>> r_to_mat_int(
    Rcpp::Matrix<Rcpp::traits::r_sexptype_traits<int>::rtype> x);

std::vector<std::vector<double>> r_to_mat_double(
    Rcpp::Matrix<Rcpp::traits::r_sexptype_traits<double>::rtype> x);

template <class T>
std::vector<std::vector<std::vector<T>>> r_to_array(Rcpp::List x)
{
    int n_elements = x.size();
    std::vector<std::vector<std::vector<T>>> x_mat(n_elements);

    for (size_t i = 0; i < n_elements; i++)
    {
        Rcpp::List x_i(x[i]);
        int nrows = x_i.size();
        x_mat[i] = std::vector<std::vector<T>>(nrows);
        for (size_t j = 0; j < nrows; j++)
        {
            Rcpp::NumericVector x_i_j(x_i[j]);
            x_mat[i][j] = Rcpp::as<std::vector<T>>(x_i_j);
        }
    }

    return x_mat;
};

std::vector<std::vector<std::vector<bool>>> r_to_array_bool(Rcpp::List x);

std::vector<std::vector<std::vector<int>>> r_to_array_int(Rcpp::List x);

std::vector<std::vector<std::vector<double>>> r_to_array_double(Rcpp::List x);

template <class T>
void print(T x)
{
    Rcpp::Rcout << x << "\n";
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
    R_FlushConsole();
#endif
};

template <class T>
void rewrite_line(T x)
{
    Rcpp::Rcout << "\r" << x;
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
    R_FlushConsole();
#endif
}

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
