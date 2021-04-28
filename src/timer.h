#ifndef TIMER_H
#define TIMER_H

#include <chrono>

template <typename TimeT = std::chrono::milliseconds,
          typename ClockT = std::chrono::steady_clock>
class Timer
{
   private:
    using timep_t = typename ClockT::time_point;
    timep_t start_ = ClockT::now();
    timep_t end_ = {};

   public:
    void tick()
    {
        end_ = timep_t{};
        start_ = ClockT::now();
    }

    void tock() { end_ = ClockT::now(); }

    template <typename TT = TimeT>
    TT duration() const
    {
        return std::chrono::duration_cast<TT>(end_ - start_);
    }
};

#endif /* TIMER_H */
