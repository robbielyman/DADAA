// PluginDADAA.hpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"
#include "cmath"

namespace DADAA {

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
  /* // second anti-derivative
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
  */
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
    // differential operator
    inline double d4(double in1, double in2, double in3, double in4, double in5) {
      double delta = in1 - in5;
      if (abs(delta) > mEps) {
        return 4.f * (d3(in1, in2, in3, in4) - d3(in2, in3, in4, in5)) / delta;
      }
      else {
        double barx = 0.5f * (in1 + in5);
      double part1 = 24.f * (waveshape4(in2) / ((barx - in2) * (barx - in2) * (in2 - in3))
        + waveshape4(in4) / ((barx - in4) * (barx - in4) * (in3 - in4))) / (in2 - in4);
      double part2 = 24.f * waveshape4(in3) / ((barx - in3) * (barx - in3) * (in2 - in3) * (in3 - in4));
      double part3 = 24.f * (waveshape3(barx) 
        - waveshape4(barx) * (1/(barx - in2) + 1/(barx - in3) + 1/(barx - in4))) 
        / ((barx - in2) * (barx - in3) * (barx - in4));
        return part1 + part2 + part3;
      }
    }
    // differential operator
    inline double d3(double in1, double in2, double in3, double in4) {
      double delta = in1 - in4;
      return 3.f * (d2(in1, in2, in3) - d2(in2, in3, in4)) / delta;
    }
    // differential operator
    inline double d2(double in1, double in2, double in3) {
      double delta = in1 - in3;
      return 2.f * (d1(in1, in2) - d1(in2, in3)) / delta;
    }
    // difference quotient
    inline double d1(double in1, double in2) {
      double delta = in1 - in2;
      return (waveshape4(in1) - waveshape4(in2)) / delta;
    }
    // Calc function
    void next_k(int nSamples);
    void next_a(int nSamples);

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
    /*// second anti-derivative
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
  */
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
    // returns true if all samples are either less than 1.f in absolute value
    // or if all samples are greater than 1.f in absolute value
    inline bool check(double nn, double n, double c, double l, double ll) {
      bool check1 = (abs(nn) <= 1.f) && (abs(n) <= 1.f) && (abs(c) <= 1.f) && (abs(l) <= 1.f) && (abs(ll) <= 1.f);
      bool check2 = (nn >= 1.f) && (n >= 1.f) && (c >= 1.f) && (l >= 1.f) && (ll >= 1.f);
      bool check3 = (nn <= -1.f) && (n <= -1.f) && (c <= -1.f) && (l <= -1.f) && (ll <= -1.f);
      return check1 || check2 || check3;
    }
    // differential operator
    inline double d4(double in1, double in2, double in3, double in4, double in5) {
      double delta = in1 - in5;
      if (abs(delta) > mEps) {
        return 4.f * (d3(in1, in2, in3, in4) - d3(in2, in3, in4, in5)) / delta;
      }
      else {
        double barx = 0.5f * (in1 + in5);
      double part1 = 24.f * (waveshape4(in2) / ((barx - in2) * (barx - in2) * (in2 - in3))
        + waveshape4(in4) / ((barx - in4) * (barx - in4) * (in3 - in4))) / (in2 - in4);
      double part2 = 24.f * waveshape4(in3) / ((barx - in3) * (barx - in3) * (in2 - in3) * (in3 - in4));
      double part3 = 24.f * (waveshape3(barx) 
        - waveshape4(barx) * (1/(barx - in2) + 1/(barx - in3) + 1/(barx - in4))) 
        / ((barx - in2) * (barx - in3) * (barx - in4));
        return part1 + part2 + part3;
      }
    }
    // differential operator
    inline double d3(double in1, double in2, double in3, double in4) {
      double delta = in1 - in4;
      return 3.f * (d2(in1, in2, in3) - d2(in2, in3, in4)) / delta;
    }
    // differential operator
    inline double d2(double in1, double in2, double in3) {
      double delta = in1 - in3;
      return 2.f * (d1(in1, in2) - d1(in2, in3)) / delta;
    }
    // difference quotient
    inline double d1(double in1, double in2) {
      double delta = in1 - in2;
      return (waveshape4(in1) - waveshape4(in2)) / delta;
    }
    // Calc function
    void next_k(int nSamples);
    void next_a(int nSamples);

    // Member variables
    const double mEps = 0.0001;
    double mNextNext = 0;
    double mNext = 0;
    double mCurrent = 0;
    double mLast = 0;
    double mLastLast = 0;
};

} // namespace DADAA
