#ifndef MCMC_PROGRESS_BAR_H
#define MCMC_PROGRESS_BAR_H

#include "timer.h"

#include <Rcpp.h>

#include <progress_bar.hpp>

class MCMCProgressBar : public ProgressBar
{
   public:
    MCMCProgressBar(int burnin_, int sample_);
    ~MCMCProgressBar(){};

    void display();
    void update(float progress);
    void end_display();

   private:
    std::string time_remaining_string_(double dur, float progress);
    std::string current_ticks_string_(float progress);
    std::string construct_ticks_display_string_(int ticks);
    void finalize_display_();

    int max_ticks_ = 50;
    int burnin_;
    int sample_;
    Timer<> clock_;
    bool finalized_ = false;
    bool timer_flag_ = false;
};

#endif /* MCMC_PROGRESS_BAR_H */
