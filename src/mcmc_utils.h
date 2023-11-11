#pragma once

#ifndef MCMC_UTILS_H_
#define MCMC_UTILS_H_

#include <Rcpp.h>
#include <algorithm>

// #include <boost/math/distributions.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/random.hpp>
#include <boost/range/algorithm.hpp>

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

long int r_to_long_int(SEXP x);

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
    std::size_t nrow = x.nrow();
    std::size_t ncol = x.ncol();
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
    std::size_t n_elements = x.size();
    std::vector<std::vector<std::vector<T>>> x_mat(n_elements);

    for (size_t i = 0; i < n_elements; i++)
    {
        Rcpp::List x_i(x[i]);
        std::size_t nrows = x_i.size();
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
    if (v.size() == 0)
    {
        Rcpp::Rcout << "[]\n";
    }
    else
    {
        Rcpp::Rcout << '[';
        for (size_t i = 0; i < v.size() - 1; i++)
        {
            Rcpp::Rcout << v[i] << ", ";
        }
        Rcpp::Rcout << v[v.size() - 1] << "]\n";
    }
}

template <class T>
inline double logit(const T x)
{
    if (x < .5)
    {
        return std::log(x) - std::log1p(-x);
    }
    else
    {
        return std::log(x / (1 - x));
    }
}

template <class T>
inline std::vector<double> logitVec(const std::vector<T> &x)
{
    std::vector<double> res;
    std::transform(x.begin(), x.end(), std::back_inserter(res),
                   UtilFunctions::logit<T>);
    return res;
}

template <class T>
inline std::pair<std::vector<double>, std::vector<double>> log_pq(
    const std::vector<T> &x)
{
    std::vector<double> logp;
    logp.reserve(x.size());
    std::vector<double> logq;
    logq.reserve(x.size());

    for (const auto el : x)
    {
        double ex = std::exp(el);
        if (el < 0)
        {
            logq.push_back(-std::log1p(ex));
            logp.push_back(logq.back() + el);
        }
        else
        {
            logp.push_back(-std::log1p(1 / ex));
            logq.push_back(logp.back() - el);
        }
    }

    return std::pair<std::vector<double>, std::vector<double>>(logp, logq);
}

template <class T>
inline std::pair<double, double> log_pq(const T x)
{
    double ex = std::exp(x);
    double logp;
    double logq;
    if (x < 0)
    {
        logq = -std::log1p(ex);
        logp = logq + x;
    }
    else
    {
        logp = -std::log1p(1 / ex);
        logq = logp - x;
    }

    return std::pair<double, double>(logp, logq);
}

inline double logitSum(const std::vector<double> &x)
{
    auto x_sorted = x;
    std::sort(x_sorted.rbegin(), x_sorted.rend());
    auto lpq = log_pq(x_sorted);

    double out;
    double cumsum = 0;

    if (x_sorted[0] < 0)
    {
        double lp1 = lpq.first[0];
        for (std::size_t i = 1; i < lpq.first.size(); ++i)
        {
            cumsum += std::exp(lpq.first[i] - lp1);
        }
        out = lp1 + std::log1p(cumsum);
    }
    else
    {
        double lq1 = lpq.second[0];
        for (std::size_t i = 1; i < lpq.first.size(); ++i)
        {
            cumsum += std::exp(lpq.first[i]);
        }
        out = std::log1p(-std::exp(lq1) + cumsum);
    }

    return out;
}

inline std::vector<double> logitScale(std::vector<double> &x, double scale)
{
    std::vector<double> out;
    out.reserve(x.size());
    std::vector<double> l2;
    l2.reserve(x.size());

    std::vector<double> u;
    u.reserve(x.size());
    std::vector<double> v;
    v.reserve(x.size());
    std::vector<double> ev;
    ev.reserve(x.size());
    std::vector<double> eumo;
    eumo.reserve(x.size());

    bool ok;
    for (std::size_t ii = 0; ii < x.size(); ++ii)
    {
        ok = (scale < std::log(2)) and
             (std::abs(scale) < std::abs(x[ii] + scale));
        u.push_back(-scale - (!ok) * x[ii]);
        v.push_back(-scale - ok * x[ii]);
        ev.push_back(std::exp(v.back()));
        eumo.push_back(std::expm1(u.back()));

        if (std::isinf(eumo.back()))
        {
            l2.push_back(std::max(u.back(), v.back()) +
                         std::log1p(std::exp(-std::abs(u.back() - v.back()))));
        }
        else
        {
            l2.push_back(std::log(eumo.back() + ev.back()));
        }

        if (v.back() > std::log(2 * std::abs(eumo.back())))
        {
            out.push_back(-(v.back() + std::log1p(eumo.back() / ev.back())));
        }
        else
        {
            out.push_back(-l2.back());
        }
    }

    return out;
}

template <class T>
inline double expit(const T x)
{
    return 1 / (1 + std::exp(-x));
}

template <class T>
inline std::vector<T> expitVec(std::vector<T> x)
{
    std::vector<double> out;
    out.reserve(x.size());
    std::transform(x.begin(), x.end(), std::back_inserter(out),
                   UtilFunctions::expit<T>);
    return out;
}

template <typename Engine>
std::vector<int> randomSequence(int min, int max, Engine rng)
{
    assert(min < max);
    std::vector<int> indices(max - min);
    std::iota(std::begin(indices), std::end(indices), min);

    auto int_generator = [=](int max_val)
    {
        boost::random::uniform_int_distribution<> dist_{0, max_val - 1};
        return dist_(*rng);
    };

    boost::range::random_shuffle(indices, int_generator);
    return indices;
}

inline double logSumExp(const double a, const double b)
{
    double max_el = std::max(a, b);
    if (max_el == -std::numeric_limits<double>::infinity())
    {
        return -std::numeric_limits<double>::infinity();
    }
    double sum = std::exp(a - max_el) + std::exp(b - max_el);
    return max_el + std::log(sum);
}

/**
 * Numerically stable log(∑(exp(a)))
 * Does not assume the iterable is sorted.
 * @tparam Iter implements iterable
 * @param begin iterator pointer to beginning
 * @param end iterator pointer to end
 * @return
 */
template <typename Iter>
typename std::iterator_traits<Iter>::value_type logSumExp(const Iter &begin,
                                                          const Iter &end)
{
    using ValueType = typename std::iterator_traits<Iter>::value_type;

    if (begin == end)
    {
        return ValueType{};
    }

    auto max_el = *std::max_element(begin, end);

    if (max_el == -std::numeric_limits<ValueType>::infinity())
    {
        return -std::numeric_limits<ValueType>::infinity();
    }

    auto sum = std::accumulate(begin, end, ValueType{},
                               [max_el](ValueType a, ValueType b)
                               { return a + std::exp(b - max_el); });
    return max_el + std::log(sum);
}

template <typename It,
          typename T = std::decay_t<decltype(*begin(std::declval<It>()))>>
double logSumExp(const It &x)
{
    return logSumExp(x.begin(), x.end());
}

template <typename T>
constexpr const T &clamp(T &el, T &low, T &high)
{
    return el < low ? low : el > high ? high : el;
}

}  // namespace UtilFunctions

#endif  // MCMC_UTILS_H_
