//
//  crtbp.cpp
//
//
//  Created by Protoss Probe on 2017/06/07.
//  Copyright © 2016-2017年 probe. All rights reserved.
//

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;
using namespace boost::math;
using namespace boost::numeric::odeint;

class CRTBP {
  public:
    CRTBP() = default;
    ~CRTBP() = default;
    CRTBP(double mu) : mu(mu) {}
    const double mu = 0.001;

  private:
};
int main() { cout << "Test Output!" << endl; }