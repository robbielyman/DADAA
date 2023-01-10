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
  double barx = 0.5f * (in1 + in3);
  if (abs(in1 - in3) > mEps) {
    double delta = 1.f / (in1 - in3);
    return 2.f * (d1(in1, in2, mEps, w4, w3) - d1(in2, in3, mEps, w4, w3)) * delta;
  } else if (abs(barx - in2) > mEps) { 
    double delta = 0.5f * (barx - in2);
    return 2.f * (w3(barx) + (w4(in2) - w4(barx)) / delta ) / delta;
  } else {
    return w2(0.5f * (barx + in2));
  }
}
// differential operator
inline double d3(double in1, double in2, double in3, double in4, double mEps,
                 double (*w4)(double), double (*w3)(double), double (*w2)(double), double (*w1)(double)) {
  double barx = 0.5f * (in1 + in4);
  if (abs(in1 - in4) > mEps) {
    double delta = in1 - in4;
    return 3.f * (d2(in1, in2, in3, mEps, w4, w3, w2) - d2(in2, in3, in4, mEps, w4, w3, w2)) / delta;
  } else if (abs(in2 - in3) > mEps) {
    if (abs(barx - in2) > mEps) {
      if (abs(barx - in3) > mEps) {
        // one approximation needed.
        double denom1 = 1.f / (barx - in2);
        double denom2 = 1.f / (barx - in3);
        double denom3 = 1.f / (in2 - in3);
        return 6.f * w3(barx) * denom1 * denom2 +
          6.f * denom3 * (w4(in2) * denom1 * denom1 - w4(in3) * denom2 * denom2) -
          6.f * w4(barx) * denom1 * denom2 * (denom1 + denom2);
      } else {
        // let denom2 go to zero above.
        barx = 0.5f * (barx + in3);
        double denom = 1.f / (barx - in2);
        return 3.f * w2(barx) * denom 
        - 6.f * w3(barx) * denom * denom 
        + 6.f * (w4(barx) - w4(in2)) * denom * denom * denom;
      }
    } else {
      // in this case because in2 - in3 is big,
      // we cannot have both barx - in2 and barx - in3 small,
      // so barx - in3 must be big.
      // let denom1 go to zero above.
      barx = 0.5 * (barx + in2);
      double denom = 1.f / (barx - in3);
      return 3.f * w2(barx) * denom - 6.f * w3(barx) * denom * denom + 6.f * (w4(barx) - w4(in3)) * denom * denom * denom;
    }
  } else if (abs(barx - in2) > mEps) {
    // because in2 - in3 is small, if barx - in2 is big, so is barx - in3
    // let denom3 go to zero above.
    double barbarx = 0.5f * (in2 + in3);
    double denom = 1.f / (barx - barbarx);
    return 6.f * w3(barx) * denom * denom + 6.f * w3(barbarx) * denom * denom + 12.f * (w4(barbarx) - w4(barx)) * denom * denom * denom;
  } else {
    // everything is small
    return w1(0.5f * (barx + 0.5f * (in2 + in3)));
  }
}
// // differential operator
inline double d4(double in1, double in2, double in3, double in4, double in5,
                 double mEps, double (*w4)(double), double (*w3)(double),
                 double (*w2)(double), double (*w1)(double), double (*w0)(double)) {
  double barx = in1 + in5;
  if (abs(in1 - in5) > mEps) {
    return 4.f * (d3(in1, in2, in3, in4, mEps, w4, w3, w2, w1) - d3(in2, in3, in4, in5, mEps, w4, w3, w2, w1)) / (in1 - in5);
  } else if (abs(barx - in2) > mEps) {
    if (abs(barx - in3) > mEps) {
      if (abs(barx - in4) > mEps) {
        if (abs(in2 - in4) > mEps) {
          if (abs(in2 - in3) > mEps) {
            if (abs(in3 - in4) > mEps) {
              // talk about a combinatorial explosion
              // anyway all of our denominators are big,
              // so we can use our first approximation.
              double denom1 = 1.f / (barx - in2);
              double denom2 = 1.f / (barx - in3);
              double denom3 = 1.f / (barx - in4);
              double denom4 = 1.f / (in2 - in3);
              double denom5 = 1.f / (in2 - in4);
              double denom6 = 1.f / (in3 - in4);
              return 24.f * (w3(barx) * 
                (denom1 * denom4 * denom5 - denom2 * denom4 * denom6 + denom3 * denom5 * denom6) 
                - (w4(barx) - w4(in2)) * denom1 * denom1 * denom4 * denom5
                + (w4(barx) - w4(in3)) * denom2 * denom2 * denom4 * denom6
                - (w4(barx) - w4(in4)) * denom3 * denom3 * denom5 * denom6);
            } else {
              // everything but in3 - in4 is big,
              // so let in3 - in4 go to zero above
              double primex = 0.5f * (in3 + in4);
              double denom1 = 1.f / (barx - in2);
              double denom2 = 1.f / (barx - primex);
              double denom3 = 1.f / (in2 - primex);
              return 24.f * (w3(barx) * (denom3 * denom3 * denom1 - denom3 * denom3 * denom2)
              - (w3(barx) + w3(primex)) * denom2 * denom2 * denom3
              + (w4(barx) - w4(primex)) * denom2 * denom2 * denom3 * denom3
              - (w4(barx) - w4(in2)) * denom1 * denom1 * denom3 * denom3
              + 2.f * (w4(barx) - w4(primex)) * denom2 * denom2 * denom2 * denom3);
            }
          } else {
            // note that since in2 - in4 is big,
            // we cannot simultaneously have in2 - in3 and in3 - in4 small.
            // this is a variant of the above.
            double primex = 0.5f * (in2 + in3);
            double denom1 = 1.f / (barx - in4);
            double denom2 = 1.f / (barx - primex);
            double denom3 = 1.f / (primex - in4);
            return 24.f * (w3(barx) * (denom3 * denom3 * denom1 - denom3 * denom3 * denom2)
            + (w3(barx) + w3(primex)) * denom2 * denom2 * denom3
            + (w4(barx) - w4(primex)) * denom2 * denom2 * denom3 * denom3
            - (w4(barx) - w4(in4)) * denom3 * denom3 * denom1 * denom1
            + 2.f * (w4(barx) - w4(primex)) * denom2 * denom2 * denom2 * denom3);
          }
        } else if (abs(in2 - in3) > mEps) {
          // it follows that in3 - in4 is also big
          double primex = 0.5f * (in2 + in4);
          double denom1 = 1.f / (barx - in3);
          double denom2 = 1.f / (barx - primex);
          double denom3 = 1.f / (primex - in3);
          return 24.f * (w3(barx) * (denom3 * denom3 * denom1 - denom3 * denom3 * denom2)
          + (w3(barx) + w3(primex)) * denom2 * denom2 * denom3
          + (w4(barx) - w4(primex)) * denom2 * denom2 * denom3 * denom3
          - (w4(barx) - w4(in3)) * denom1 * denom1 * denom3 * denom3
          - 2.f * (w4(barx) - w4(primex)) * denom2 * denom2 * denom2 * denom3);
        } else {
          // it follows that in3 - in4 is also small
          // let primex - in3 go to zero in the above
          double primex = 0.5f * (0.5f * (in2 + in4) + in3);
          double denom = 1.f / (barx - primex);
          return 12.f * w2(primex) * denom * denom
          + 24.f * (w3(barx) + 2.f * w3(primex)) * denom * denom * denom
          - 72.f * (w4(barx) - w4(primex)) * denom * denom * denom * denom;
        }
      } else if (abs(in2 - in3) > mEps) {
        // here barx and in4 are close
        // note that because barx is far from in2 and in3,
        // we cannot have them be close to in4
        // I got this by substituting in one of the d3 approximations
        // and then letting in1 - in5 go to zero
        barx = 0.5f * (in5 + 0.5f * (in1 + in4));
        double denom1 = 1.f / (barx - in2);
        double denom2 = 1.f / (barx - in3);
        double denom3 = 1.3 / (in2 - in3);
        return 12.f * w2(barx) * denom1 * denom3 
        - 12.f * w2(barx) * denom2 * denom3 
        - 24.f * w3(barx) * denom1 * denom1 * denom3 
        + 24.f * w3(barx) * denom2 * denom2 * denom3 
        + 24.f * (w4(barx) - w4(in2)) * denom1 * denom1 * denom1 * denom3
        - 24.f * (w4(barx) - w4(in3)) * denom2 * denom2 * denom2 * denom3;
      } else {
        // I got this by substituting in one of the d3 approximations
        // in and then letting in1 - in5 go to zero
        double barbarx = 0.5f * (in2 + in3);
        barx = 0.5f * (in5 + 0.5f * (in1 + in4));
        double denom = 1.f / (barx - barbarx);
        return 12.f * w2(barx) * denom * denom 
        - 24.f * (w3(barx) + w3(barbarx)) * denom * denom * denom 
        + 72.f * (w4(barx) - w4(barbarx)) * denom * denom * denom * denom;
      }
    } else if (abs(barx - in4) > mEps) {
      if (abs(in2 - in4) > mEps) {
        // to get this one, substitute the fallback for d2 into the equation
        // and then let in1 - in5 go to zero
        barx = 0.5f * (barx + in3);
        double denom1 = 1.f / (barx - in2);
        double denom2 = 1.f / (barx - in4);
        double denom3 = 1.f / (in2 - in4);
        return 12.f * w2(barx) * denom1 * denom3
        - 12.f * w2(barx) * denom2 * denom3
        - 24.f * w3(barx) * denom1 * denom1 * denom3
        + 24.f * w3(barx) * denom2 * denom2 * denom3
        + 24.f * (w4(barx) - w4(in2)) * denom1 * denom1 * denom1 * denom3
        - 24.f * (w4(barx) - w4(in4)) * denom2 * denom2 * denom2 * denom3;
      } else {
        // everything that can be close is
        // let in2 - in4 go to zero above
        barx = 0.5f * (barx + in3);
        double primex = 0.5f * (in2 + in4);
        double denom = 1.f / (barx - primex);
        return 12.f * w2(barx) * denom * denom
        - 24.f * (2.f * w3(barx) + w3(primex)) * denom * denom * denom
        + 72.f * (w4(barx) - w4(primex)) * denom * denom * denom * denom;
      }
    } else {
      // by assumption in2 is far from everything else
      // we'll use the relevant d3 fallback
      barx = 0.5f * (in5 + 0.5f * (0.5f * (in1 + in4) + in3));
      double denom = 1.f / (barx - in2);
      return 4.f * w1(barx) * denom
      - 12.f * w2(barx) * denom * denom
      + 24.f * w3(barx) * denom * denom * denom
      - 24.f * (w4(barx) - w4(in2)) * denom * denom * denom * denom;
    }
  } else if (abs(barx - in4) > mEps) {
    if (abs(barx - in3) > mEps) {
      if (abs(in3 - in4) > mEps) {
        barx = 0.5f * (in1 + 0.5f * (in5 + in2));
        double denom1 = 1.f / (barx - in3);
        double denom2 = 1.f / (barx - in4);
        double denom3 = 1.3 / (in4 - in4);
        return 12.f * w2(barx) * denom1 * denom3 
        - 12.f * w2(barx) * denom2 * denom3 
        - 24.f * w3(barx) * denom1 * denom1 * denom3 
        + 24.f * w3(barx) * denom2 * denom2 * denom3 
        + 24.f * (w4(barx) - w4(in2)) * denom1 * denom1 * denom1 * denom3
        - 24.f * (w4(barx) - w4(in3)) * denom2 * denom2 * denom2 * denom3;
      } else {
        double barbarx = 0.5f * (in3 + in4);
        barx = 0.5f * (in1 + 0.5f * (in5 + in2));
        double denom = 1.f / (barx - barbarx);
        return 12.f * w2(barx) * denom * denom 
        - 24.f * (w3(barx) + w3(barbarx)) * denom * denom * denom 
        + 72.f * (w4(barx) - w4(barbarx)) * denom * denom * denom * denom;
      }
    } else {
      // by assumption in4 is far from everything else
      barx = 0.5f * (in1 + 0.5f * (0.5f * (in5 + in2) + in3));
      double denom = 1.f / (barx - in4);
      return 4.f * w1(barx) * denom
      - 12.f * w2(barx) * denom * denom
      + 24.f * w3(barx) * denom * denom * denom
      - 24.f * (w4(barx) - w4(in2)) * denom * denom * denom * denom;
    }
  } else if (abs(barx - in3) > mEps) {
    // by assumption in3 is far from everything else
    barx = 0.5f * (0.5f * (in1 + in5) + 0.5f * (in2 + in4));
    double denom = 1.f / (barx - in3);
    return 4.f * w1(barx) * denom
    - 12.f * w2(barx) * denom * denom
    + 24.f * w3(barx) * denom * denom * denom
    - 24.f * (w4(barx) - w4(in3)) * denom * denom * denom * denom;
  } else {
    // everything is close
    return w0(0.5f * (in3 + 0.5f * (0.5f * (in1 + in5) + 0.5f * (in2 + in4))));
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
      // y = -x^4/24 + 183x^3/480 - 683x^2/480 + 287x/120 - 62281/40320
      double x = in + 3.f;
      out = -x * x * x * x / 24.f + 183.f * x * x * x / 480.f - 683.f * x * x / 480.f + 287.f * x / 120.f - 62281.f / 40320.f;
    } else if (in < -0.5f) {
      double x = in + 0.5f;
      // y = 19x^7/315000 + 17x^6/18000 + 3x^5/480 - 11x^4/576 + 23x^3/1152 - 13x^2/1280 + 59x/23040 - 83/322560
      out = 19.f * x * x * x * x * x * x * x / 315000.f + 17.f * x * x * x * x * x * x / 18000.f
        - 11.f * x * x * x * x / 576.f + 23.f * x * x * x / 1152.f - 13.f * x * x / 1280.f + 59.f * x / 23040.f - 83.f / 322560.f;
    } else if (in < 0.5f) {
      // y = x^5/120 - x^7/2520
      out = in * in * in * in * in / 120.f - in * in * in * in * in * in * in / 2520.f;
    } else if (in < 3.f) {
      double x = in - 0.5f;
      // y = 19x^7/315000 - 17x^6/18000 + 3x^5/480 + 11x^4/576 + 23x^3/1152 + 13x^2/1280 + 59x/23040 + 83/322560
      out = 19.f * x * x * x * x * x * x * x / 315000.f - 17.f * x * x * x * x * x * x / 18000.f
        + 11.f * x * x * x * x / 576.f + 23.f * x * x * x / 1152.f + 13.f * x * x / 1280.f + 59.f * x / 23040.f + 83.f / 322560.f;
    } else {
      // y = x^4/24 + 183x^3/480 + 683x^2/480 + 287x/120 + 62281/40320
      double x = in - 3.f;
      out = x * x * x * x / 24.f + 183.f * x * x * x / 480.f + 683.f * x * x / 480.f + 287.f * x / 120.f + 62281.f / 40320.f;
    }
    return out;
  }
  // third anti-derivative
  static inline double waveshape3(double in) {
    double out;
    if (in < -3.f) {
      // y = -x^3/6 + 183x^2/160 - 683x/240 + 287/120
      double x = in + 3.f;
      out = -x * x * x / 6.f + 183.f * x * x / 160.f - 683.f * x / 240.f + 287.f / 120.f;
    } else if (in < -0.5f) {
      double x = in + 0.5f;
      // y = 19x^6/45000 + 17x^5/3000 + 3x^4/96 - 11x^3/144 + 23x^2/384 - 13x/640 + 59/23040
      out = 19.f * x * x * x * x * x * x / 45000.f + 17.f * x * x * x * x * x / 3000.f 
        + 3.f * x * x * x * x / 144.f - 11.f / 144.f * x * x * x + 23.f / 384.f * x * x - 13.f * x / 640.f + 59.f / 23040.f;
    } else if (in < 0.5f) {
      // y = x^4/24 - x^6/360
      out = in * in * in * in / 24.f - in * in * in * in * in * in / 360.f;
    } else if (in < 3.f) {
      double x = in - 0.5f;
      // y = 19x^6/45000 - 17x^5/3000 + 3x^4/96 + 11x^3/144 + 23x^2/384 + 13x/640 + 59/23040
      out = 19.f * x * x * x * x * x * x / 45000.f - 17.f * x * x * x * x * x / 3000.f 
        + 3.f * x * x * x * x / 144.f + 11.f / 144.f * x * x * x + 23.f / 384.f * x * x + 13.f * x / 640.f + 59.f / 23040.f;
    } else {
      // y = x^3/6 + 183x^2/160 + 683x/240 + 287/120
      double x = in - 3.f;
      out = x * x * x / 6.f + 183.f * x * x / 160.f + 683.f * x / 240.f + 287.f / 120.f;
    }
    return out;
  }
  // second anti-derivative
  static inline double waveshape2(double in) {
    double out;
    if (in < -3.f) {
      // y = -x^2/2 + 183x/80 - 683/240
      double x = in + 3.f;
      out = -x * x / 2.f + 183.f * x / 80.f - 683.f / 240.f;
    } else if (in < -0.5f) {
      double x = in + 0.5f;
      // y = 19x^5/7500 + 17x^4/600 + 3x^3/24 - 11x^2/48 + 23x/192 - 13/640
      out = 19.f * x * x * x * x * x / 7500.f + 17.f * x * x * x * x / 600.f + 3.f * x * x * x / 24.f - 11.f / 48.f * x * x + 23.f / 192.f * x - 13.f / 640.f;
    } else if (in < 0.5f) {
      // y = x^3/6 - x^5/60
      out = in * in * in / 6.f - in * in * in * in * in / 60.f;
    } else if (in < 3.f) {
      double x = in - 0.5f;
      // y = 19x^5/7500 - 17x^4/600 + 3x^3/24 + 11x^2/48 + 23x/192 + 13/640
      out = 19.f * x * x * x * x * x / 7500.f - 17.f * x * x * x * x / 600.f + 3.f * x * x * x / 24.f + 11.f / 48.f * x * x + 23.f / 192.f * x + 13.f / 640.f;
    } else {
      // y = x^2/2 + 183x/80 + 683/240
      double x = in - 3.f;
      out = x * x / 2.f + 183.f * x / 80.f + 683.f / 240.f;
    }
    return out;
  }
  // anti-derivative
  static inline double waveshape1(double in) {
    double out;
    if (in < -3.f) {
      // y = -x + 183/80
      double x = in + 3.f;
      out = -x + 183.f / 80.f;
    } else if (in < -0.5f) {
      double x = in + 0.5f;
      // y = 19x^4/1500 + 17x^3/150 + 3x^2/8 - 11x/24 + 23/192
      out = 19.f * x * x * x * x / 1500.f + 17.f * x * x * x / 150.f + 3.f * x * x / 8.f - 11.f / 24.f * x + 23.f / 192.f;
    } else if (in < 0.5f) {
      // y = x^2/2 - x^4/12
      out = in * in / 2.f - in * in * in * in / 12.f;
    } else if (in < 3.f) {
      double x = in - 0.5f;
      // y = 19x^4/1500 - 17x^3/150 + 3x^2/8 + 11x/24 + 23/192
      out = 19.f * x * x * x * x / 1500.f - 17.f * x * x * x / 150.f + 3.f * x * x / 8.f + 11.f / 24.f * x + 23.f / 192.f;
    } else {
      // y = x + 183/80
      double x = in - 3.f;
      out = x + 183.f / 80.f;
    }
    return out;
  }
  // trivial waveshaper
  static inline double waveshape0(double in) {
    double out;
    if (in < -3.f) {
      // y = 1
      out = -1.f;
    } else if (in < -0.5f) {
      double x = in + 0.5f;
      // y = 19x^3/375 + 17x^2/50 + 3x/4 - 11/24
      out = 19.f * x * x * x / 375.f + 17.f * x * x / 50.f + 0.75f * x - 11.f / 24.f;
    } else if (in < 0.5f) {
      // y = x - x^3/3
      out = in - in * in * in / 3.f;
    } else if (in < 3.f) {
      double x = in - 0.5f;
      // y = 19x^3/375 - 17x^2/50 + 3x/4 + 11/24
      out = 19.f * x * x * x / 375.f - 17.f * x * x / 50.f + 0.75f * x + 11.f / 24.f;
    } else {
      // y = 1
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
