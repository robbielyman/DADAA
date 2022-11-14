// PluginDADAA.cpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#include "SC_PlugIn.hpp"
#include "DADAA.hpp"

static InterfaceTable* ft;

namespace DADAA {

CADAA::CADAA() {
  if (isAudioRateIn(1)) {
    mCalcFunc = make_calc_function<CADAA, &CADAA::next_a>();
    next_a(1);
  } else {
    mCalcFunc = make_calc_function<CADAA, &CADAA::next_k>();
    next_k(1);
  }
}

TADAA::TADAA() {
  if (isAudioRateIn(1)) {
    mCalcFunc = make_calc_function<TADAA, &TADAA::next_a>();
    next_a(1);
  } else {
    mCalcFunc = make_calc_function<TADAA, &TADAA::next_k>();
    next_k(1);
  }
}

void CADAA::next_k(int nSamples) {
    const float* input = in(0);
    const float gain = in0(1);
    float* wet = out(0);
    float* dry = out(1);

    // anti-aliased hard clipping function
    for (int i = 0; i < nSamples; ++i) {
      mLastLast = mLast;
      mLast = mCurrent;
      mCurrent = mNext;
      mNext = mNextNext;
      mNextNext = input[i];
      dry[i] = mCurrent;
      // if all five samples are in one "regime", just hard-clip
      if (check(mNextNext, mNext, mCurrent, mLast, mLastLast)) {
        wet[i] = waveshape0(mCurrent * gain);
      }
      // otherwise do ADAA process
      else {
        wet[i] = d4(mNextNext * gain, mNext * gain, mCurrent * gain, mLast * gain, mLastLast * gain);
      }
    }
}

void CADAA::next_a(int nSamples) {
    const float* input = in(0);
    const float* gain = in(1);
    float* wet = out(0);
    float* dry = out(1);

    // anti-aliased hard clipping function
    for (int i = 0; i < nSamples; ++i) {
      mLastLast = mLast;
      mLast = mCurrent;
      mCurrent = mNext;
      mNext = mNextNext;
      mNextNext = input[i];
      dry[i] = mCurrent;
      // if all five samples are in one "regime", just hard-clip
      if (check(mNextNext, mNext, mCurrent, mLast, mLastLast)) {
        wet[i] = waveshape0(mCurrent * gain[i]);
      }
      // otherwise do ADAA process
      else {
        wet[i] = d4(mNextNext * gain[i], mNext * gain[i], mCurrent * gain[i], mLast * gain[i], mLastLast * gain[i]);
      }
    }
}

void TADAA::next_k(int nSamples) {
    const float* input = in(0);
    const float gain = in0(1);
    float* wet = out(0);
    float* dry = out(1);

    // anti-aliased hard clipping function
    for (int i = 0; i < nSamples; ++i) {
      mLastLast = mLast;
      mLast = mCurrent;
      mCurrent = mNext;
      mNext = mNextNext;
      mNextNext = input[i];
      dry[i] = mCurrent;
      wet[i] = d4(mNextNext * gain, mNext * gain, mCurrent * gain, mLast * gain, mLastLast * gain);
    }
}

void TADAA::next_a(int nSamples) {
    const float* input = in(0);
    const float* gain = in(1);
    float* wet = out(0);
    float* dry = out(1);

    // anti-aliased hard clipping function
    for (int i = 0; i < nSamples; ++i) {
      mLastLast = mLast;
      mLast = mCurrent;
      mCurrent = mNext;
      mNext = mNextNext;
      mNextNext = input[i];
      dry[i] = mCurrent;
      wet[i] = d4(mNextNext * gain[i], mNext * gain[i], mCurrent * gain[i], mLast * gain[i], mLastLast * gain[i]);
    }
}

} // namespace DADAA

PluginLoad(DADAAUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<DADAA::CADAA>(ft, "CADAA", false);
    registerUnit<DADAA::TADAA>(ft, "TADAA", false);
}
