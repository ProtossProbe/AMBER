//
//  crtbp.cpp
//
//
//  Created by Protoss Probe on 2017/06/07.
//  Copyright © 2016-2017年 probe. All rights reserved.
//

#include "crtbp.hpp"
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace boost::math;
using namespace boost::numeric::odeint;

int main() {
    double mu = 0.001;
    pcrtbp system(mu, 'p');
    double a_res = pow(0.5, 2. / 3.);
    system.print_vec(system.elements2vector({{1, 0, 0, 0}}));
}