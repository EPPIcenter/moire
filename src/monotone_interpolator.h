#pragma once

#ifndef MONOTONE_INTERPOLATOR_H_
#define MONOTONE_INTERPOLATOR_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

// Shape-preserving cubic Hermite interpolation through knots with strictly
// increasing x and monotone y (Fritsch & Carlson 1980, with the
// Fritsch-Butland slope formula used by PCHIP). If the knot values are
// monotone the interpolant is monotone. Evaluation is clamped to the knot
// range. MCMC::adapt_temp() uses this to invert the cumulative communication
// barrier: temperature as a function of barrier.
class MonotoneCubicInterpolator
{
   public:
    // Requires x.size() == y.size() >= 2 and x strictly increasing.
    MonotoneCubicInterpolator(std::vector<double> x, std::vector<double> y)
        : x_(std::move(x)), y_(std::move(y)), m_(x_.size(), 0.0)
    {
        const std::size_t n = x_.size();
        std::vector<double> h(n - 1);
        std::vector<double> delta(n - 1);
        for (std::size_t k = 0; k + 1 < n; ++k)
        {
            h[k] = x_[k + 1] - x_[k];
            delta[k] = (y_[k + 1] - y_[k]) / h[k];
        }

        if (n == 2)
        {
            m_[0] = m_[1] = delta[0];
            return;
        }

        for (std::size_t k = 1; k + 1 < n; ++k)
        {
            if (delta[k - 1] * delta[k] <= 0.0)
            {
                m_[k] = 0.0;  // local extremum or flat: keep it flat
            }
            else
            {
                const double w1 = 2.0 * h[k] + h[k - 1];
                const double w2 = h[k] + 2.0 * h[k - 1];
                m_[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k]);
            }
        }
        m_[0] = edge_slope(h[0], h[1], delta[0], delta[1]);
        m_[n - 1] = edge_slope(h[n - 2], h[n - 3], delta[n - 2], delta[n - 3]);
    }

    double operator()(double x) const
    {
        if (x <= x_.front())
        {
            return y_.front();
        }
        if (x >= x_.back())
        {
            return y_.back();
        }
        const std::size_t k =
            static_cast<std::size_t>(std::upper_bound(x_.begin(), x_.end(), x) - x_.begin()) - 1;
        const double h = x_[k + 1] - x_[k];
        const double t = (x - x_[k]) / h;
        const double t2 = t * t;
        const double t3 = t2 * t;
        const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
        const double h10 = t3 - 2.0 * t2 + t;
        const double h01 = -2.0 * t3 + 3.0 * t2;
        const double h11 = t3 - t2;
        return h00 * y_[k] + h10 * h * m_[k] + h01 * y_[k + 1] + h11 * h * m_[k + 1];
    }

   private:
    // Three-point one-sided estimate for an end knot, clamped so the end
    // segment stays shape preserving. h0/d0 belong to the end segment.
    static double edge_slope(double h0, double h1, double d0, double d1)
    {
        double m = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
        if (m * d0 <= 0.0)
        {
            return 0.0;
        }
        if (d0 * d1 <= 0.0 && std::fabs(m) > 3.0 * std::fabs(d0))
        {
            return 3.0 * d0;
        }
        return m;
    }

    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> m_;
};

#endif  // MONOTONE_INTERPOLATOR_H_
