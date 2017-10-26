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
#include <ctime>
#include <fstream>
#include <iostream>
#include <libiomp/omp.h>
#include <string>

using namespace std;
using namespace boost::math;
using namespace boost::numeric::odeint;

int main(int argc, char *argv[]) {
    string inputstring(argv[1]);
    string infostring(argv[2]);
    auto input = readInputFromTxt(inputstring);
    auto info = readInfoFromTxt(infostring);
    int number = input.size();
    cout << "Total Number: " << number << endl;

    double endt = info[0], dt = info[1];
    int jump = info[2];

    crtbp system;
    orbit3d orbits[MAX_NUMBER];
    for (size_t i = 0; i < number; i++) {
        orbits[i].setElement(input[i].first);
        orbits[i].setState(crtbp::elementsToRot(orbits[i].getElement(), 0));
        orbits[i].setDt(dt);
        orbits[i].setName(to_string(i + 1));
    }
    double start = omp_get_wtime();
    system.inteNbody(orbits, number, endt, jump);
    cout << "Time: " << omp_get_wtime() - start << endl;
}