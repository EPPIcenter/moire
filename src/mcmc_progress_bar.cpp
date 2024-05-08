#include "mcmc_progress_bar.h"

#include "mcmc_utils.h"

MCMCProgressBar::MCMCProgressBar(int burnin, int sample)
    : burnin_(burnin), sample_(sample)
{
}

void MCMCProgressBar::display()
{
    using namespace UtilFunctions;
    print("0%    10   20   30   40   50   60   70   80   90   100%");
    print("[----|----|----|----|----|----|----|----|----|----]");
};

void MCMCProgressBar::update(float progress)
{
    if (finalized_) return;

    if (timer_flag_)
    {
        timer_flag_ = false;
        clock_.tick();
    }
    else
    {
        // stop and record time no more than every .1 seconds
        if (clock_.time_since_event(events::UPDATE_CONSOLE).count() < 1000)
        {
            return;
        }

        clock_.tock();
        auto duration = clock_.duration();

        std::string time_string =
            time_remaining_string_(duration.count(), progress);
        std::string ticks_string = current_ticks_string_(progress);
        std::string llik_string = current_llik_string_();
        std::string empty_space = std::string(time_string.length(), ' ');

        std::stringstream ss;
        ss << "|" << ticks_string << "| " << time_string << llik_string
           << empty_space;

        // clear console line and update
        clock_.record_event(events::UPDATE_CONSOLE);

        #ifdef __EMSCRIPTEN__
        UtilFunctions::print(ss.str());
        #else
        UtilFunctions::rewrite_line(ss.str());
        #endif

        if (progress == 1)
        {
            finalize_display_();
        }
    }
};

void MCMCProgressBar::end_display() { update(1); };

std::string MCMCProgressBar::time_remaining_string_(float dur, float progress)
{
    float rem_time = (dur / progress) * (1 - progress);
    int hour = 0;
    int min = 0;
    int sec = 0;

    hour = rem_time / (60 * 60 * 1000);
    rem_time = std::fmod(rem_time, (60 * 60 * 1000));

    min = rem_time / (60 * 1000);
    rem_time = std::fmod(rem_time, (60 * 1000));

    sec = rem_time / 1000;

    std::stringstream ss;
    if (hour != 0) ss << hour << "h ";
    if (min != 0 or hour != 0) ss << min << "m ";
    ss << sec << "s ";
    std::string time_str = ss.str();
    return time_str;
};

std::string MCMCProgressBar::current_ticks_string_(float progress)
{
    int ticks_remaining = (int)(progress * max_ticks_);
    return construct_ticks_display_string_(ticks_remaining);
}

std::string MCMCProgressBar::construct_ticks_display_string_(int ticks)
{
    std::stringstream ss;

    int burnin_tick =
        (int)(((float)burnin_ / (burnin_ + sample_)) * max_ticks_);

    for (int i = 0; i < max_ticks_ - 1; ++i)
    {
        if (i == burnin_tick)
        {
            ss << "B";
        }
        else if (i < ticks)
        {
            ss << "*";
        }
        else
        {
            ss << " ";
        }
    }
    return ss.str();
}

std::string MCMCProgressBar::current_llik_string_()
{
    std::stringstream ss;
    ss << "(Llik: " << llik_ << ", Hot Chain: " << hot_chain_ << ")";
    return ss.str();
}

void MCMCProgressBar::finalize_display_()
{
    if (finalized_) return;

    UtilFunctions::print("");
    finalized_ = true;
}

void MCMCProgressBar::set_llik(float llik) { llik_ = llik; }
void MCMCProgressBar::set_hot_chain(int idx) { hot_chain_ = idx; }
