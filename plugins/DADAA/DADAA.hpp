// PluginDADAA.hpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"
#include "cmath"

namespace DADAA {

// difference quotient
inline double d1(double in1, double in2, double mEps, double (*w4)(double), double (*w3)(double)) {
  double delta = in1 - in2;
  if (abs(delta) > mEps) {
    return (w4(in1) - w4(in2)) / delta;
  } else {
    return w3(0.5f * (in1 + in2));
  }
}
// differential operator
inline double d2(double in1, double in2, double in3, double mEps,
                 double (*w4)(double), double (*w3)(double), double (*w2)(double)) {
  double delta = in1 - in3;
  double barx = 0.5f * (in1 + in3);
  double deldelta = barx - in2;
  if (abs(delta) > mEps) {
    return 2.f * (d1(in1, in2, mEps, w4, w3) - d1(in2, in3, mEps, w4, w3)) / delta;
  } else if (abs(deldelta) > mEps) {
    return 2.f * (w3(barx) + (w4(in2) - w4(barx)) / deldelta) / deldelta;
  } else {
    return w2(0.5f * (barx + in2));
  }
}
// differential operator
inline double d3(double in1, double in2, double in3, double in4, double mEps,
                 double (*w4)(double), double (*w3)(double), double (*w2)(double), double (*w1)(double)) {
  double delta = in1 - in4;
  double deldelta = in2 - in3;
  double barx = 0.5f * (in1 + in4);
  double barbarx = 0.5f * (in2 + in3);
  double deldeldelta = barbarx - barx;
  if (abs(delta) > mEps) {
    return 3.f * (d2(in1, in2, in3, mEps, w4, w3, w2) - d2(in2, in3, in4, mEps, w4, w3, w2)) / delta;
  } else if (abs(deldelta) > mEps) {
    double part1 = 6.f * ((w4(in2) - w4(barx)) / ((in2 - barx) * (in2 - barx)) - (w4(in3) - w4(barx)) / ((in3 - barx) * (in3 - barx))) / deldelta;
    double part2 = 6.f * w3(barx) / ((in2 - barx) * (in3 - barx));
    return part1 + part2;
  } else if (abs(deldeldelta) > mEps) {
    return 6.f * (w3(barx) + w3(barbarx)) / (deldeldelta * deldeldelta) + 12.f * (w4(barbarx) - w4(barx)) / (deldeldelta * deldeldelta * deldeldelta);
  } else {
    return w1(0.5f * (barx + barbarx));
  }
}
// differential operator
inline double d4(double in1, double in2, double in3, double in4, double in5, 
                 double mEps, double (*w4)(double), double (*w3)(double), 
                 double (*w2)(double), double (*w1)(double), double (*w0)(double)) {
  double delta = in1 - in5;
  double deldelta = in2 - in4;
  double barx = 0.5f * (in1 + in5);
  double barbarx = 0.5f * (in2 + in4);
  double deldeldelta = barbarx - barx;
  double primex = 0.5f * (barbarx + barx);
  double deldeldeldelta = primex - in3;
  if (abs(delta) > mEps) {
    return 4.f * (d3(in1, in2, in3, in4, mEps, w4, w3, w2, w1) 
      - d3(in2, in3, in4, in5, mEps, w4, w3, w2, w1)) / delta;
  } else if (abs(deldelta) > mEps) {
    double part1 = 24.f * (w4(in2) / ((barx - in2) * (barx - in2) * (in2 - in3))
      + w4(in4) / ((barx - in4) * (barx - in4) * (in3 - in4))) / (in2 - in4);
    double part2 = 24.f * w4(in3) / ((barx - in3) * (barx - in3) * (in2 - in3) * (in3 - in4));
    double part3 = 24.f * (w3(barx) - w4(barx) * (1/(barx - in2) + 1/(barx - in3) + 1/(barx - in4)))
      / ((barx - in2) * (barx - in3) * (barx - in4));
    return part1 - part2 + part3;
  } else if (abs(deldeldelta) > mEps) {
    double part1 = 24.f * (w3(barbarx) / (barbarx - in3) + w3(barx) / (barx - in3) 
      + (w4(barbarx) - w4(barx)) / ((barx - in3) * (barx - in3)) 
      - 2.f * (w4(barbarx) - w4(barx)) / ((barx - in3) * (barbarx - barx)));
    double part2 = 24.f * (w4(in3) - w4(barbarx)) / ((barbarx - in3) * (barbarx - in3) * (barx - in3) * (barx - in3));
    return part1 + part2;
  } else if (abs(deldeldeldelta) > mEps) {
    double part1 = 4.f * w1(primex) / deldeldeldelta;
    double part2 = 12.f * w2(primex) / (deldeldeldelta * deldeldeldelta);
    double part3 = 24.f * w3(primex) / (deldeldeldelta * deldeldeldelta * deldeldeldelta);
    double part4 = 24.f * (w4(in3) - w4(primex)) / (deldeldeldelta * deldeldeldelta * deldeldeldelta * deldeldeldelta);
    return part1 - part2 + part3 + part4;
  } else {
    return w0(0.5f * (primex + in3));
  }
}

class TADAA : public SCUnit {
public:
  TADAA();
    
  // Destructor
  // ~TADAA();
private:
  // final anti-derivative
  static inline double waveshape4(double in) {
    double out;
    if (in < -3.f) {
      double x = in + 3.f;
      out = -35.f * x * x * x * x + 980.f * x * x * x / 3.f - 4781.f * x * x / 4.f
        + 2009.f * x - 62281.f / 48.f;
    } else if (in < -0.5f) {
      double x = in + 0.5f;
      out = 19.f * x * x * x * x * x * x * x / 375.f + 119.f * x * x * x * x * x * x / 150.f
        + 21.f * x * x * x * x * x / 4.f - 385.f * x * x * x * x / 24.f 
        + 805.f * x * x * x / 48.f - 273.f * x * x / 32.f + 413.f * x / 192.f - 83.f / 384.f;
    } else if (in < 0.5f) {
      out = 7.f * in * in * in * in * in - in * in * in * in * in * in *in / 3.f;
    } else if (in < 3.f) {
      double x = in - 0.5f;
      out = 19.f * x * x * x * x * x * x * x / 375.f - 119.f * x * x * x * x * x * x / 150.f
        + 21.f * x * x * x * x * x / 4.f + 385.f * x * x * x * x / 24.f
        + 805.f * x * x * x / 48.f + 273.f * x * x / 32.f + 413.f * x / 192.f + 83.f / 384.f;
    } else {
      double x = in - 3.f;
      out = 35.f * x * x * x * x + 980.f * x * x * x / 3.f + 4781.f * x * x / 4.f 
        + 2009.f * x + 62281.f / 48.f;
    }
    return out / 840.f;
  }
  // third anti-derivative
  static inline double waveshape3(double in) {
    double out;
    if (in < -3.f) {
      double x = in + 3.f;
      out = -20.f * x * x * x + 140.f * x * x - 683.f * x / 2.f + 287.f;
    } else if (in < -0.5f) {
      double x = in + 0.5f;
      out = 19.f * x * x * x * x * x * x / 375.f + 17.f * x * x * x * x * x / 25.f
        + 15.f * x * x * x * x / 4.f - 55.f * x * x * x / 6.f + 115.f * x * x / 16.f
        - 39.f * x / 16.f + 59.f / 192.f;
    } else if (in < 0.5f) {
      out = 5.f * in * in * in * in - in * in * in * in * in * in / 3.f;
    } else if (in < 3.f) {
      double x = in - 0.5f;
      out = 19.f * x * x * x * x * x * x / 375.f - 17.f * x * x * x * x * x / 25.f
        + 15.f * x * x * x * x / 4.f + 55.f * x * x * x / 6.f + 115.f * x * x / 16.f
        + 39.f * x / 16.f + 59.f / 192.f;
    } else {
      double x = in - 3.f;
      out = 20.f * x * x * x + 140.f * x * x + 683.f * x / 2.f + 287.f;
    }
    return out / 120.f;
  }
  // second anti-derivative
  static inline double waveshape2(double in) {
    double out;
    if (in < -3.f) {
      double x = in + 3.f;
      out = -60.f * x * x + 280.f * x - 683.f / 2.f;
    } else if (in < -0.5f) {
      double x = in + 0.5f;
      out = 38.f * x * x * x * x * x / 125.f + 17.f * x * x * x * x / 5.f + 15.f * x * x * x - 55.f * x * x / 2.f + 115.f * x / 8.f - 39.f / 16.f;
    } else if (in < 0.5f) {
      out = 20.f * in * in * in - 2.f * in * in * in * in * in;
    } else if (in < 3.f) {
      double x = in - 0.5f;
      out = 38.f * x * x * x * x * x / 125.f - 17.f * x * x * x * x / 5.f + 15.f * x * x * x + 55.f * x * x / 2.f + 115.f * x / 8.f + 39.f / 16.f;
    } else {
      double x = in - 3.f;
      out = 60.f * x * x + 280.f * x + 683.f / 2.f;
    }
    return out / 120.f;
  }
  // anti-derivative
  static inline double waveshape1(double in) {
    double out;
    if (in < -3.f) {
      double x = in + 3.f;
      return -24.f * x + 56.f;
    } else if (in < -0.5f) {
      double x = in + 0.5f;
      out = 38.f * x * x * x * x / 125.f + 68.f * x * x * x / 25.f + 9.f * x * x - 11.f * x + 23.f / 8.f;
    } else if (in < 0.5f) {
      out = 12.f * in * in - 2.f * in * in * in * in;
    } else if (in < 3.f) {
      double x = in - 0.5f;
      out = 38.f * x * x * x * x / 125.f - 68.f * x * x * x / 25.f + 9.f * x * x + 11.f * x + 23.f / 8.f;
    } else {
      double x = in - 3.f;
      out = 24.f * x + 56.f;
    }
    return out / 24.f;
  }
  // trivial waveshaper
  static inline double waveshape0(double in) {
    double out;
    if (in < -3.f) {
      out = -1.f;
    } else if (in < -0.5f) {
      double x = in + 0.5f;
      out = 19.f * x * x * x / 375.f + 17.f * x * x / 50.f + 0.75f * x - 11.f / 24.f;
    } else if (in < 0.5f) {
      out = in - in * in * in / 3.f;
    } else if (in < 3.f) {
      double x = in - 0.5f;
      out = 19.f * x * x * x / 375.f + 17.f * x * x / 50.f + 0.75f * x + 11.f / 24.f;
    } else {
      out = 1.f;
    }
    return out;
  }
    // Calc function
    void next(int nSamples);

    // Member variables
    const double mEps = 0.0001;
    double mNextNext = 0;
    double mNext = 0;
    double mCurrent = 0;
    double mLast = 0;
    double mLastLast = 0;
};

class CADAA : public SCUnit {
public:
    CADAA();

    // Destructor
    // ~CADAA();

private:
    // final anti-derivative
    static inline double waveshape4(double in) {
      double out;
      if (in < -1.f) {
        double x = in + 1.f;
        out = -5.f * x * x * x * x;
      }
      else if (in < 1.f) {
        out = in * in * in * in * in - 10.f * in * in * in - 20.f * in * in - 15.f * in - 4.f;
      }
      else {
        double x = in - 1.f;
        out = 5.f * x * x * x * x - 40.f * x * x - 80.f * x - 48.f;
      }
      return out / 120.f;
    }
    // third anti-derivative
    static inline double waveshape3(double in) {
      double out;
      if (in < -1.f) {
        double x = in + 1.f;
        out = -4.f * x * x * x;
      }
      else if (in < 1.f) {
        out = in * in * in * in - 6.f * in * in - 8.f * in - 3.f;
      }
      else {
        double x = in - 1.f;
        out = 4.f * x * x * x - 16.f * x - 16;
      }
      return out / 24.f;
    }
    // second anti-derivative
    static inline double waveshape2(double in) {
      double out;
      if (in < -1.f) {
        double x = in + 1.f;
        out = -3.f * x * x;
      }
      else if (in < 1.f) {
        out = in * in * in - 3.f * in - 2.f;
      }
      else {
        double x = in - 1.f;
        out = 3.f * x * x - 4.f;
      }
      return out / 6.f;
    }
    // first anti-derivative
    static inline double waveshape1(double in) {
      double out;
      if (in < -1.f) {
        double x = in + 1.f;
        out = -2.f * x;
      }
      else if (in < 1.f) {
        out = in * in - 1.f;
      }
      else {
        double x = in - 1.f;
        out = 2.f * x;
      }
      return out / 2.f;
    }
    // trivial waveshaper
    static inline double waveshape0(double in) {
      if (in < -1.f) {
        return -1.f;
      }
      else if (in < 1.f) {
        return in;
      }
      else {
        return 1.f;
      }
    }
    // Calc function
    void next(int nSamples);

    // Member variables
    const double mEps = 0.00001;
    double mNextNext = 0;
    double mNext = 0;
    double mCurrent = 0;
    double mLast = 0;
    double mLastLast = 0;
};

} // namespace DADAA
