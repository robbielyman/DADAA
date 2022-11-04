// PluginCADAA.hpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"
#include "cmath"

namespace CADAA {

class CADAA : public SCUnit {
public:
    CADAA();

    // Destructor
    // ~CADAA();

private:
    // final anti-derivative
    inline double waveshape4(double in) {
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
    inline double waveshape3(double in) {
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
    inline double waveshape2(double in) {
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
    inline double waveshape1(double in) {
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
    inline double waveshape0(double in) {
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
      bool check2 = (abs(nn) >= 1.f) && (abs(n) >= 1.f) && (abs(c) >= 1.f) && (abs(l) >= 1.f) && (abs(ll) >= 1.f);
      return check1 || check2;
    }
    // difference quotient
    inline double dMinus(double in1, double in2) {
      double delta = in1 - in2;
      // if not ill-conditioned, form difference quotient
      if (delta > mEps) {
        return (waveshape4(in1) - waveshape4(in2)) / delta;
      }
      // otherwise, use approximation of difference quotient
      else {
        double avg = 0.5 * (in1 + in2);
        return waveshape3(avg);
      }
    }
    // difference quotient, used in approximation
    inline double dMinusShift(double in1, double in2) {
      double delta = in1 - in2;
      // if not ill-conditioned, form difference quotient
      if (delta > mEps) {
        return (waveshape3(in1) - waveshape3(in2)) / delta;
      }
      // otherwise, use approximation of difference quotient
      else {
        double avg = 0.5 * (in1 + in2);
        return waveshape2(avg);
      }
    }
    // difference quotient, used in approximation
    inline double dMinusShiftShift(double in1, double in2) {
      double delta = in1 - in2;
      // if not ill-conditioned, form difference quotient
      if (delta > mEps) {
        return (waveshape2(in1) - waveshape2(in2)) / delta;
      }
      // otherwise, use approximation of difference quotient
      else {
        double avg = 0.5 * (in1 + in2);
        return waveshape1(avg);
      }
    }
    // differential operator
    inline double d2(double in1, double in2, double in3) {
      double delta = in1 - in3;
      // if not ill-conditioned, apply differential operator
      if (delta > mEps) {
        return 2.f * (dMinus(in1, in2) - dMinus(in2, in3)) / delta;
      }
      // otherwise, use approximation of differential operator
      else {
        double avg = 0.5 * (in1 + in3);
        double deldelta = avg - in2;
        if (deldelta > mEps) {
          return 2.f * (waveshape3(avg) + (waveshape4(in2) - waveshape4(avg)) / deldelta) / deldelta;
        }
        else {
          double avavg = 0.5 * (avg + in2);
          return waveshape2(avavg);
        }
      }
    }
    // differential operator, used in approximation
    inline double d2Shift(double in1, double in2, double in3) {
      double delta = in1 - in3;
      // if not ill-conditioned, apply differential operator
      if (delta > mEps) {
        return 2.f * (dMinusShift(in1, in2) - dMinusShift(in2, in3)) / delta;
      }
      // otherwise, use approximation of differential operator
      else {
        double avg = 0.5 * (in1 + in3);
        double deldelta = avg - in2;
        if (deldelta > mEps) {
          return 2.f * (waveshape2(avg) + (waveshape3(in2) - waveshape3(avg)) / deldelta) / deldelta;
        }
        else {
          double avavg = 0.5 * (avg + in2);
          return waveshape1(avavg);
        }
      }
    }
    // differential operator, used in approximation
    inline double d2ShiftShift(double in1, double in2, double in3) {
      double delta = in1 - in3;
      // if not ill-conditioned, apply differential operator
      if (delta > mEps) {
        return 2.f * (dMinusShiftShift(in1, in2) - dMinusShiftShift(in2, in3)) / delta;
      }
      // otherwise, use approximation of differential operator
      else {
        double avg = 0.5 * (in1 + in3);
        double deldelta = avg - in2;
        if (deldelta > mEps) {
          return 2.f * (waveshape1(avg) + (waveshape2(in2) - waveshape2(avg)) / deldelta) / deldelta;
        }
        else {
          double avavg = 0.5 * (avg + in2);
          return waveshape0(avavg);
        }
      }
    }
    // differential operator
    inline double dMinusD2(double in1, double in2, double in3, double in4) {
      double delta = in2 - in3;
      // if not ill-conditioned, apply differential operator
      if (delta > mEps) {
        return (d2(in1, in2, in3) - d2(in2, in3, in4)) / delta;
      }
      // otherwise, use approximation of differential operator
      else {
        double av1 = 0.5 * (in1 + in2);
        double av2 = 0.5 * (in2 + in3);
        double av3 = 0.5 * (in3 + in4);
        return d2Shift(av1, av2, av3);
      }
    }
    // ADAA process
    inline double d2D2(double in1, double in2, double in3, double in4, double in5) {
      double delta = in2 - in4;
      // if not ill-conditioned, apply differential operator
      if (delta > mEps) {
        return 2.f * (dMinusD2(in1, in2, in3, in4) - dMinusD2(in2, in3, in4, in5)) / delta;
      }
      // otherwise, use approximation of differential operator
      else {
        double av1 = 0.5 * (in1 + in3);
        double av2 = 0.5 * (in2 + in4);
        double av3 = 0.5 * (in3 + in5);
        double deldelta = av2 - in3;
        if (deldelta > mEps) {
          return 2.f * (d2Shift(av1, av2, av3) + (d2(in2, in3, in4) - d2(av1, av2, av3)) / deldelta) / deldelta;
        }
        else {
          double avav1 = 0.5 * (av1 + in2);
          double avav2 = 0.5 * (av2 + in3);
          double avav3 = 0.5 * (av3 + in4);
          return d2ShiftShift(avav1, avav2, avav3);
        }
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

} // namespace CADAA
