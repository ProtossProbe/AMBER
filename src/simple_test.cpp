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
#include <string>
#include <sstream>

using namespace std;
using namespace boost::math;
using namespace boost::numeric::odeint;

int main()
{
    ofstream outputfile;
    outputfile.open("../assets/output.txt");
    state_type4 x = {{1, 2, 3, 4}};
    outputfile << setprecision(12) << x[0] << '\t' << x[1] << endl;
    outputfile.close();
}