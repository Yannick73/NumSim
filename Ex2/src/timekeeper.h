#pragma once

#include <cassert>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <mpi.h>

inline std::chrono::time_point<std::chrono::steady_clock> timestamp()
{
    std::chrono::time_point<std::chrono::steady_clock> t{std::chrono::steady_clock::now()};
    return t;
}

inline double getDurationS(std::chrono::time_point<std::chrono::steady_clock> t0)
{
    std::chrono::time_point<std::chrono::steady_clock> tEnd = timestamp();
    assert(tEnd > t0);
    const auto elapsed_time{tEnd - t0};
    return ((double)elapsed_time.count()) * 1e-9;
}

//! Small time class, probably overengineered for the job
class Timekeeper
{
public:
    //! updates the communication time, time needs to be bigger or equal 0.0 and is measured in s
    inline void addTime(double time) { assert(time >= 0.0); cummulatedTime_ += time; }

    //! gets the cummulated communication time
    inline double getCummulatedTime() const { return cummulatedTime_; }

    //! sets the last time point, which can be used
    inline void setT0() { t0_ = timestamp(); }

    //! gets the duration since t0 in seconds
    inline double getDurationSinceT0s() { return getDurationS(t0_); }

    //! adds the time since t0 to the current cummulated time
    inline void addTimeSinceT0() { addTime(getDurationSinceT0s()); }

    //! gets the last set t0
    inline std::chrono::time_point<std::chrono::steady_clock> getT0() { return t0_; }

private:
    //! last t0
    std::chrono::time_point<std::chrono::steady_clock> t0_ = timestamp();

    //! all time spent tracked by the time-keeper
    double cummulatedTime_ = 0.0;
};
