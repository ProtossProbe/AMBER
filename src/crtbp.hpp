//
//  crtbp.cpp
//
//
//  Created by Protoss Probe on 2017/06/07.
//  Copyright © 2016-2017年 probe. All rights reserved.
//

#ifndef _CRTBP_HPP_
#define _CRTBP_HPP_

#include "crtbp.hpp"
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>
typedef boost::array<double, 2> vec2;
typedef boost::array<double, 4> vec4;
typedef boost::array<double, 6> vec6;
typedef boost::array<double, 31> vec31;
typedef boost::array<double, 61> vec61;

const double pi = M_PI;

class pcrtbp {
  public:
    pcrtbp();
    ~pcrtbp();
    pcrtbp(double mu);
    pcrtbp(double mu, char option);
    const double mu = 0.001;
    const char option = 'p';
    vec4 vector = {{0, 0, 0, 0}};
    vec4 elements = {{0, 0, 0, 0}};
    std::ofstream outputfile1, outputfile2, outputfile3, outputfile4,
        outputfile5;

    class planar_crtbp_ode {
      public:
        double mu;
        planar_crtbp_ode(double mu);
        void operator()(const vec4 &x, vec4 &dxdt, double t);
    };
    double jacobi_inte(const vec4 &x);
    double jacobi_inte_ele(const vec4 &elements);
    void integ_output(vec4 x, const double endtime, const double dt = 0.01);
    void write_keypoint(const vec4 &x, const double t);
    void write_detail(const vec4 &x, const double t);
    vec4 rot2init(vec4 x, const double t);
    vec4 init2rot(vec4 x, const double t);
    double rot2r(vec4 x);
    vec4 vector2elements(const vec4 &x);
    double true2mean(double theta, double e);
    vec4 elements2vector(const vec4 &elements);
    void print_vec(const vec4 &x);
    void cal_init_vec();
    double F1_func(const vec4 &x, const double t);
    double average(double a, double e, int limit);
    void test();
    double delta21(double a);
    double Hamiltion21(double delta, double Phi, double phi);
    double deltax1x2(double x1, double x2);
};

class crtbp {
  public:
    crtbp();
    ~crtbp();
};

#endif