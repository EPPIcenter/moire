#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <map>
#include <string>

template <typename eventT, typename TimeT = std::chrono::milliseconds,
          typename ClockT = std::chrono::steady_clock>
class Timer
{
   private:
    using timep_t = typename ClockT::time_point;
    timep_t start_ = ClockT::now();
    timep_t end_ = {};
    std::map<eventT, timep_t> event_tracker_{};

   public:
    void tick()
    {
        end_ = timep_t{};
        start_ = ClockT::now();
    }

    void tock() { end_ = ClockT::now(); }

    void record_event(const eventT event)
    {
        event_tracker_[event] = ClockT::now();
    }

    template <typename TT = TimeT>
    TT duration() const
    {
        return std::chrono::duration_cast<TT>(end_ - start_);
    }

    template <typename TT = TimeT>
    TT time_since_event(const eventT event)
    {
        return std::chrono::duration_cast<TT>(ClockT::now() -
                                              event_tracker_[event]);
    }
};

#endif /* TIMER_H */
