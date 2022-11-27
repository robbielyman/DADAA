// PluginDADAA.cpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#include "SC_PlugIn.hpp"
#include "DADAA.hpp"

static InterfaceTable* ft;

namespace DADAA {

CADAA::CADAA() {
  mCalcFunc = make_calc_function<CADAA, &CADAA::next>();
  next(1);
}

TADAA::TADAA() {
  mCalcFunc = make_calc_function<TADAA, &TADAA::next>();
  next(1);
}

void CADAA::next(int nSamples) {
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
    wet[i] = DADAA::d4(mNextNext * gain, mNext * gain, mCurrent * gain, mLast * gain, mLastLast * gain, 
                       mEps, waveshape4, waveshape3, waveshape2, waveshape1, waveshape0);
  }
}

void TADAA::next(int nSamples) {
    const float* input = in(0);
    const float gain = in0(1);
    float* wet = out(0);
    float* dry = out(1);

    // anti-aliased tanh saturation function
    for (int i = 0; i < nSamples; ++i) {
      mLastLast = mLast;
      mLast = mCurrent;
      mCurrent = mNext;
      mNext = mNextNext;
      mNextNext = input[i];
      dry[i] = mCurrent;
      wet[i] = DADAA::d4(mNextNext * gain, mNext * gain, mCurrent * gain, mLast * gain, mLastLast * gain,
                         mEps, waveshape4, waveshape3, waveshape2, waveshape1, waveshape0);
    }
}

} // namespace DADAA

PluginLoad(DADAAUGens) {
  // Plugin magic
  ft = inTable;
  registerUnit<DADAA::CADAA>(ft, "CADAA", false);
  registerUnit<DADAA::TADAA>(ft, "TADAA", false);
}
