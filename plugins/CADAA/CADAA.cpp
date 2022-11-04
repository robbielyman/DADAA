// PluginCADAA.cpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#include "SC_PlugIn.hpp"
#include "CADAA.hpp"

static InterfaceTable* ft;

namespace CADAA {

CADAA::CADAA() {
    mCalcFunc = make_calc_function<CADAA, &CADAA::next>();
    next(1);
}

void CADAA::next(int nSamples) {
    const float* input = in(0);
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
        wet[i] = waveshape0(mCurrent);
      }
      // otherwise do ADAA process
      else {
        wet[i] = d2D2(mNextNext, mNext, mCurrent, mLast, mLastLast);
      }
    }
}


} // namespace CADAA

PluginLoad(CADAAUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<CADAA::CADAA>(ft, "CADAA", false);
}
