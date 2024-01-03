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
    void set_llik(float llik);
    void set_hot_chain(int idx);

   private:
    enum events
    {
        UPDATE_CONSOLE
    };

    std::string time_remaining_string_(float dur, float progress);
    std::string current_ticks_string_(float progress);
    std::string current_llik_string_();
    std::string construct_ticks_display_string_(int ticks);
    void finalize_display_();

    int max_ticks_ = 50;
    int burnin_;
    int sample_;
    float llik_;
    int hot_chain_;
    Timer<events> clock_;
    bool finalized_ = false;
    bool timer_flag_ = false;
};

#endif /* MCMC_PROGRESS_BAR_H */
