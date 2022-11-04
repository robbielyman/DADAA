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
    const float* gain = in(1);
    float* outbuf = out(0);

    // simple gain function
    for (int i = 0; i < nSamples; ++i) {
        outbuf[i] = input[i] * gain[i];
    }
}

} // namespace CADAA

PluginLoad(CADAAUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<CADAA::CADAA>(ft, "CADAA", false);
}
