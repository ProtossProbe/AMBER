//
//  crtbp.cpp
//
//
//  Created by Protoss Probe on 2017/06/07.
//  Copyright © 2016-2017 probe. All rights reserved.
//

#ifndef _AVERAGE_HPP_
#define _AVERAGE_HPP_

#include "crtbp.hpp"

class resonance
{
public:
  resonance() = default;
  ~resonance() = default;
  resonance(bool isPrograde, size_t k, size_t kj, double N, double Sz,
            double Min, double Max, size_t numY, size_t num)
      : isPrograde(isPrograde), k(k), kj(kj), N(N), Sz(Sz), Min(Min),
        Max(Max), numY(numY), num(num){};
  resonance(bool isPrograde, size_t k, size_t kj, double N, double a, int phi,
            double Max, size_t numY, size_t num)
      : isPrograde(isPrograde), k(k), kj(kj), N(N), a(a), phi((double)phi),
        phi_c((double)phi), Max(Max), numY(numY), num(num){};
  resonance(double iMax, double a, size_t numY, size_t num)
      : iMax(iMax), a(a), numY(numY), num(num){};

  char option = 'm';
  double a = 1, e = 0, I = 180;
  double N = -2, Sz = 0, iMax = 0;
  double Min = 0, Max = 1, S = 0;
  double H0 = 0;
  double ome = 0, Ome = 0, M = 0, Mj = 0, phi = 0;
  double f = 0;
  const bool isPrograde = false;
  const double k = 1, kj = 1;
  const double ord = isPrograde ? abs(k - kj) : (k + kj);
  const double phi_c = 0;
  const size_t numY = 360, num = 360;

  void calAEI();
  void calAction();
  void calH0();
  void calPhi();
  void calM();
  void calTrue();
  void displayInfo();

  static double dotDisturb(const vec3 &v, const vec3 &r1);
  double calDotDisturb();
  static double ringDisturb(const vec3 &v,
                            const double ap);
  double calRingDisturb(std::string planets);

  double inteOneCycle();
  double inteTwoCycle(double phiAmp);
  double inteRing(std::string planets);

  void averagePhi();
  void averageOme(double phiAmp);
  void averageOmeRing(std::string planets);
};

#endif