#pragma once

#ifndef MCMC_UTILS_H_
#define MCMC_UTILS_H_

#include <Rcpp.h>
#include <algorithm>
#include <span>

// #include <boost/math/distributions.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/random.hpp>
#include <boost/range/algorithm.hpp>

#include "multivector.h"
#include "env_defs.h"
#include "genotyping_data.h"
#include "sampler.h"


#define OVERFLO 1e100
#define UNDERFLO 1e-100

namespace UtilFunctions
{
float fastlog2(float x);

float fastlog(float x);

int r_to_bool(SEXP x);

int r_to_int(SEXP x);

long int r_to_long_int(SEXP x);

float r_to_float(SEXP x);

std::string r_to_string(SEXP x);

std::vector<bool> r_to_vector_bool(SEXP x);

std::vector<int> r_to_vector_int(SEXP x);

std::vector<float> r_to_vector_float(SEXP x);

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

std::vector<std::vector<float>> r_to_mat_float(
    Rcpp::Matrix<Rcpp::traits::r_sexptype_traits<float>::rtype> x);

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

std::vector<std::vector<std::vector<float>>> r_to_array_float(Rcpp::List x);

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

template <class T>
void message(T x)
{   
    SEXP x_ = Rcpp::CharacterVector(x);
    Rcpp::message(x_);
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
inline float logit(const T x)
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
inline std::vector<float> logitVec(const std::vector<T> &x)
{
    std::vector<float> res;
    std::transform(x.begin(), x.end(), std::back_inserter(res),
                   UtilFunctions::logit<T>);
    return res;
}

template <typename InputIterator>
inline std::vector<float> logitVec(InputIterator begin, InputIterator end)
{
    std::vector<float> res;
    std::transform(begin, end, std::back_inserter(res),
                   UtilFunctions::logit<typename std::iterator_traits<InputIterator>::value_type>);
    return res;
}




template <class T>
inline std::pair<std::vector<float>, std::vector<float>> log_pq(
    const std::vector<T> &x)
{
    std::vector<float> logp;
    logp.reserve(x.size());
    std::vector<float> logq;
    logq.reserve(x.size());

    for (const auto el : x)
    {
        float ex = std::exp(el);
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

    return std::pair<std::vector<float>, std::vector<float>>(logp, logq);
}

template <class InputIterator>
inline std::pair<std::vector<float>, std::vector<float>> log_pq(InputIterator begin, InputIterator end)
{
    std::vector<float> logp;
    logp.reserve(std::distance(begin, end));
    std::vector<float> logq;
    logq.reserve(std::distance(begin, end));

    for (const auto el : std::span<typename std::iterator_traits<InputIterator>::value_type const>(begin, end))
    {
        float ex = std::exp(el);
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

    return std::pair<std::vector<float>, std::vector<float>>(logp, logq);
}

template <class T>
inline std::pair<float, float> log_pq(const T x)
{
    float ex = std::exp(x);
    float logp;
    float logq;
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

    return std::pair<float, float>(logp, logq);
}

inline float logitSum(const std::vector<float> &x)
{
    auto x_sorted = x;
    std::sort(x_sorted.rbegin(), x_sorted.rend());
    auto lpq = log_pq(x_sorted);

    float out;
    float cumsum = 0;

    if (x_sorted[0] < 0)
    {
        float lp1 = lpq.first[0];
        for (std::size_t i = 1; i < lpq.first.size(); ++i)
        {
            cumsum += std::exp(lpq.first[i] - lp1);
        }
        out = lp1 + std::log1p(cumsum);
    }
    else
    {
        float lq1 = lpq.second[0];
        for (std::size_t i = 1; i < lpq.first.size(); ++i)
        {
            cumsum += std::exp(lpq.first[i]);
        }
        out = std::log1p(-std::exp(lq1) + cumsum);
    }

    return out;
}

inline std::vector<float> logitScale(std::vector<float> &x, float scale)
{
    std::vector<float> out;
    out.reserve(x.size());
    std::vector<float> l2;
    l2.reserve(x.size());

    std::vector<float> u;
    u.reserve(x.size());
    std::vector<float> v;
    v.reserve(x.size());
    std::vector<float> ev;
    ev.reserve(x.size());
    std::vector<float> eumo;
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
inline float expit(const T x)
{
    return 1 / (1 + std::exp(-x));
}

template <class T>
inline std::vector<T> expitVec(std::vector<T> x)
{
    std::vector<float> out;
    out.reserve(x.size());
    std::transform(x.begin(), x.end(), std::back_inserter(out),
                   UtilFunctions::expit<T>);
    return out;
}


inline float logSumExp(const float a, const float b)
{
    float max_el = std::max(a, b);
    if (max_el == -std::numeric_limits<float>::infinity())
    {
        return -std::numeric_limits<float>::infinity();
    }
    float sum = std::exp(a - max_el) + std::exp(b - max_el);
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

#ifdef HAS_EXECUTION
    auto sum = std::reduce(std::execution::unseq, begin, end, ValueType{},
                           [max_el](ValueType a, ValueType b)
                           { return a + std::exp(b - max_el); });
#else
    auto sum = std::accumulate(begin, end, ValueType{},
                               [max_el](ValueType a, ValueType b)
                               { return a + std::exp(b - max_el); });
#endif
    return max_el + std::log(sum);
}

template <typename It,
          typename T = std::decay_t<decltype(*begin(std::declval<It>()))>>
float logSumExp(const It &x)
{
    return logSumExp(x.begin(), x.end());
}

inline float logSumExp(const std::vector<float> &x)
{
    float max_el = *std::max_element(x.begin(), x.end());
    if (max_el == -std::numeric_limits<float>::infinity())
    {
        return -std::numeric_limits<float>::infinity();
    }

    float sum = 0;
#pragma omp simd reduction(+ : sum)
    for (size_t i = 0; i < x.size(); ++i)
    {
        sum += std::exp(x[i] - max_el);
    }
    return max_el + std::log(sum);
}

template <typename T>
constexpr const T &clamp(const T &el, const T &low, const T &high)
{
    return el < low ? low : el > high ? high : el;
}

template <typename T>
T jaccard_similarity(std::span<int const> x, std::span<int const> y) {
    std::vector<int> intersection;
    std::vector<int> union_;
    std::set_intersection(x.begin(), x.end(), y.begin(), y.end(), std::back_inserter(intersection));
    std::set_union(x.begin(), x.end(), y.begin(), y.end(), std::back_inserter(union_));

    if (union_.size() == 0) {
        return 0.0;
    }

    T similarity = static_cast<T>(intersection.size()) / union_.size();
    
    // Log extreme values for debugging
    if (similarity < 0.0 || similarity > 1.0) {
        UtilFunctions::print("Warning: Jaccard similarity out of range:", similarity, "intersection:", intersection.size(), "union:", union_.size());
    }
    
    return similarity;
}

template <typename T>
MultiVector<T, 2> calculate_pairwise_jaccard_similarity(GenotypingData &genotyping_data) {
    MultiVector<T, 2> dist({genotyping_data.num_samples, genotyping_data.num_samples});
    dist.fill(0.0);

    for (size_t i = 0; i < genotyping_data.num_samples; ++i) {
        for (size_t j = i + 1; j < genotyping_data.num_samples; ++j) {
            T total_similarity = 0.0;
            for (size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                T locus_similarity = jaccard_similarity<T>(genotyping_data.get_observed_alleles(i, locus_idx), genotyping_data.get_observed_alleles(j, locus_idx));
                total_similarity += locus_similarity;
            }
            dist.at({i, j}) = total_similarity / genotyping_data.num_loci;
            dist.at({j, i}) = dist.at({i, j});
        }
    }
    
    return dist;
}

template <class T>
RaggedMultiVector<T, 3> calculate_clustered_allele_frequencies(GenotypingData &genotyping_data, int num_populations, Sampler &sampler)
{
    
    RaggedMultiVector<T, 3> p;
    p.resize({num_populations, genotyping_data.num_loci}, genotyping_data.num_alleles);

    MultiVector<T, 2> jaccard_similarity_matrix = genotyping_data.jaccard_similarity_matrix;

    // Calculate the number of samples to be used for each cluster
    int max_samples_per_cluster = std::ceil(static_cast<float>(genotyping_data.num_samples) / num_populations);
    std::vector<int> samples_per_population(num_populations, max_samples_per_cluster);
    samples_per_population[num_populations - 1] = genotyping_data.num_samples - std::accumulate(samples_per_population.begin(), samples_per_population.end() - 1, 0);
    

    // Track which samples have been assigned to clusters
    std::vector<bool> sample_assigned(genotyping_data.num_samples, false);
    std::vector<int> available_samples;
    available_samples.reserve(genotyping_data.num_samples);
    for (int i = 0; i < genotyping_data.num_samples; ++i) {
        available_samples.push_back(i);
    }

    for (size_t pop_idx = 0; pop_idx < num_populations; ++pop_idx) {
        
        if (available_samples.size() < samples_per_population[pop_idx]) {
            samples_per_population[pop_idx] = available_samples.size();
        }
        
        // select one random sample from the available samples
        std::vector<int> random_indices = sampler.sample_random_sequence(0, available_samples.size());
        int sample_idx = available_samples[random_indices[0]];

        // select the samples_per_population[pop_idx] samples closest to the selected sample based on the jaccard distance
        std::vector<int> closest_samples(samples_per_population[pop_idx]);
        closest_samples[0] = sample_idx;
        
        // Create a vector of distances from the selected sample to all other available samples
        std::vector<std::pair<float, int>> distances;
        for (size_t i = 0; i < available_samples.size(); ++i) {
            int other_sample = available_samples[i];
            if (other_sample != sample_idx) {
                distances.push_back({jaccard_similarity_matrix.at({sample_idx, other_sample}), other_sample});
            }
        }
        
        // Sort by distance in descending order and select the farthest ones
        std::sort(distances.begin(), distances.end(), std::greater<std::pair<float, int>>());
        for (size_t i = 1; i < samples_per_population[pop_idx]; ++i) {
            closest_samples[i] = distances[i-1].second;
        }
        

        // Remove the selected samples from available samples
        for (int selected_sample : closest_samples) {
            sample_assigned[selected_sample] = true;
            auto it = std::find(available_samples.begin(), available_samples.end(), selected_sample);
            if (it != available_samples.end()) {
                available_samples.erase(it);
            }
        }
        

        // use the closest samples to calculate the allele frequencies for the population
        for (size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
            for (size_t allele_idx = 0; allele_idx < genotyping_data.num_alleles[locus_idx]; ++allele_idx) {
                p.at({pop_idx, locus_idx, allele_idx}) = 2.0f;
            }
        }

        std::vector<int> total_alleles(genotyping_data.num_loci, 0);
        for (size_t sample_idx : closest_samples) {
            for (size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                for (size_t allele_idx = 0; allele_idx < genotyping_data.num_alleles[locus_idx]; ++allele_idx) {
                    p.at({pop_idx, locus_idx, allele_idx}) += genotyping_data.get_observed_alleles(sample_idx, locus_idx)[allele_idx]; 
                    total_alleles[locus_idx] += genotyping_data.get_observed_alleles(sample_idx, locus_idx)[allele_idx];
                }
            }
        }
        for (size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
            for (size_t allele_idx = 0; allele_idx < genotyping_data.num_alleles[locus_idx]; ++allele_idx) {
                p.at({pop_idx, locus_idx, allele_idx}) /= total_alleles[locus_idx];
            }
        }

        
    }

    return p;
}

}  // namespace UtilFunctions

#endif  // MCMC_UTILS_H_
