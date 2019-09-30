/*
1euro_filter: Adapted from http://cristal.univ-lille.fr/~casiez/1euro/OneEuroFilter.cc by Nicolas Roussel

Copyright (C) 2019 Xuning Yang

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include <chrono>
#include <cmath>
#include <memory>
#include <iostream>
#include <stdexcept>

using Clock = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::time_point<Clock>;

class LowPassFilter {
public:
  LowPassFilter(double alpha, double initval = 0.0) {
    y = s = initval;
    setAlpha(alpha);
    initialized  =  false;
  }
  ~LowPassFilter() {}

  double filter(double value) {
    double result;
    if (initialized)
      result = a * value + (1.0 - a) * s;
    else {
      result = value;
      initialized = true;
    }
    y = value;
    s = result;
    return result;
  }

  double filterWithAlpha(double value, double alpha) {
    setAlpha(alpha);
    return filter(value);
  }

  bool hasLastRawValue(void) {
    return initialized;
  }

  double lastRawValue(void) {
    return y;
  }

private:
  double y, a, s;
  bool initialized;

  void setAlpha(double alpha) {
    if (alpha <= 0.0 || alpha > 1.0)
      throw std::range_error("alpha should be in (0.0., 1.0]");
    a  =  alpha;
  }

};

// OneEuroFilter
class OneEuroFilter {
public:
  OneEuroFilter(double freq,
		double mincutoff = 1.0, double beta_ = 0.0, double dcutoff = 1.0) {
    setFrequency(freq);
    setMinCutoff(mincutoff);
    setBeta(beta_);
    setDerivateCutoff(dcutoff);
    x = std::make_unique<LowPassFilter>(alpha(mincutoff));
    dx = std::make_unique<LowPassFilter>(alpha(dcutoff));
    first_time_ = true;
  }
  ~OneEuroFilter() {}

  double filter(double value, TimePoint timestamp = TimePoint()) {
    auto current_time = timestamp;

    if (first_time_ != true) {
      auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(current_time - prev_time);

      // update the sampling frequency based on timestamps
      if (duration.count() != 0)
        freq  =  1.0 / duration.count();
    } else {
      first_time_ = false;
    }
    prev_time = current_time;

    // estimate the current variation per second
    double dvalue = x->hasLastRawValue() ? (value - x->lastRawValue()) * freq : 0.0; // FIXME: 0.0 or value?
    double edvalue = dx->filterWithAlpha(dvalue, alpha(dcutoff));

    // use it to update the cutoff frequency
    double cutoff = mincutoff + beta_ * std::abs(edvalue);

    // filter the given value
    return x->filterWithAlpha(value, alpha(cutoff));
  }

private:
  double freq;
  double mincutoff;
  double beta_;
  double dcutoff;
  std::unique_ptr<LowPassFilter> x;
  std::unique_ptr<LowPassFilter> dx;
  std::chrono::high_resolution_clock::time_point prev_time;
  bool first_time_;

  double alpha(double cutoff) {
    double te  =  1.0 / freq;
    double tau  =  1.0 / (2 * M_PI * cutoff);
    return 1.0 / (1.0 + tau/te);
  }

  void setFrequency(double f) {
    if (f <= 0) throw std::range_error("freq should be >  0");
    freq  =  f;
  }

  void setMinCutoff(double mc) {
    if (mc <= 0) throw std::range_error("mincutoff should be > 0");
    mincutoff  =  mc;
  }

  void setBeta(double b) {
    beta_  =  b;
  }

  void setDerivateCutoff(double dc) {
    if (dc <= 0) throw std::range_error("dcutoff should be > 0");
    dcutoff  =  dc;
  }
};
