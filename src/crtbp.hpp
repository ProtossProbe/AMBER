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
#include <boost/foreach.hpp>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <iostream>
#include <libiomp/omp.h>
#include <math.h>
#include <string>
#include <utility>
#include <vector>

const double pi = M_PI;
const double pi180 = pi / 180;
const double mu = 0.001;
const size_t MAX_NUMBER = 10000;
const std::string GLOBAL_OUTPUT_LOCATION = "assets/_output/";
const std::string GLOBAL_LOCATION = "assets/";
typedef boost::array<double, 2> vec2;
typedef boost::array<double, 3> vec3;
typedef boost::array<double, 4> vec4;
typedef boost::array<double, 5> vec5;
typedef boost::array<double, 6> vec6;
typedef boost::array<double, 7> vec7;
typedef boost::array<double, 12> vec12;

static vec6 joinVector3(const vec3 &v1, const vec3 &v2) {
    return {{v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]}};
}
static vec12 joinVector6(const vec6 &v1, const vec6 &v2) {
    return {{v1[0], v1[1], v1[2], v1[3], v1[4], v1[5], v2[0], v2[1], v2[2],
             v2[3], v2[4], v2[5]}};
}

static double vector3Norm(const vec3 &v) {
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static double vector6Norm(const vec6 &v) {
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3] +
                v[4] * v[4] + v[5] * v[5]);
}

static double vec6Dot(const vec6 &v1, const vec6 &v2) {
    double result = 0;
    for (size_t i = 0; i < v1.size(); i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

static vec3 vec3Cross(const vec3 &v1, const vec3 &v2) {
    vec3 v3;
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return v3;
}

class orbit3d {
  private:
    double div = 1e-8;
    double megno_temp = 0;
    double megno_fast = 0;
    double megno_sum = 0;

  public:
    std::string name = "default";
    double time = 0;
    double dt = 0.001;
    double lcn = 0;
    double megno = 0;
    double megno_max = 0;
    size_t steps = 0;
    std::ofstream outputfile;
    vec12 vec = {{0, 0, 0, 0, 0, 0, div, div, div, div, div, div}};
    vec6 ele = {{1, 0, 0, 0, 0, 0}};
    const double getTime() { return time; }
    const vec6 getState() {
        return {{vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]}};
    }

    const vec6 getDelta() {
        return {{vec[6], vec[7], vec[8], vec[9], vec[10], vec[11]}};
    }
    const vec3 getPosition() { return {{vec[0], vec[1], vec[2]}}; }
    const vec3 getVelocity() { return {{vec[3], vec[4], vec[5]}}; }
    std::string getName() { return name; };
    const vec6 getElement() { return ele; };

    void setState(const vec6 &state);
    void setElement(vec6 element) { ele = element; };
    void setName(std::string newname) { name = newname; };
    void setDt(double t) { dt = t; };
    void setOutputFile() {
        std::string location = GLOBAL_OUTPUT_LOCATION + "Ast_" + getName();
        outputfile.open(location + ".txt");
    };
    void closeOutputFile() { outputfile.close(); };
    const double getJacobi();
    vec2 deltaNormCal();

    double getLCN();
    void updateMEGNO(const double ticktime);
    double getMEGNO();
    double getMEGNOMax();
};

static void printInteData(orbit3d &orb) {
    std::cout << "State:  ";
    double j;
    BOOST_FOREACH (j, orb.getPosition())
        std::cout << std::setw(12) << j << " ";
    std::cout << "  |   "
              << "Time: " << std::setw(6) << orb.getTime()
              << "   |   Jacobi: " << std::setw(8) << orb.getJacobi()
              << "   |   LCN: " << std::setw(8) << orb.getLCN()
              << "   |   MEGNO: " << std::setw(8) << orb.getMEGNO()
              << std::endl;
}

static void writeInteData(orbit3d &orb) {
    orb.outputfile << std::setw(8) << orb.getTime() << '\t';
    double j;
    BOOST_FOREACH (j, orb.getPosition())
        orb.outputfile << std::setw(10) << j << '\t';
    orb.outputfile << std::setw(10) << orb.getJacobi() << '\t' << std::setw(10)
                   << orb.getMEGNO() << std::endl;
}

struct observer {
    orbit3d *orb;
    size_t jump;
    double ticktime;
    observer(orbit3d &orbit, size_t jump, double tick)
        : orb(&orbit), jump(jump), ticktime(tick) {}
    void operator()(const vec12 &x, double t) {
        if (t > 0) {
            (*orb).steps++;
            (*orb).time = t;
            (*orb).updateMEGNO(ticktime);
            auto megno = (*orb).getMEGNOMax();
            if (megno >= 8) {
                throw 8.0;
                return;
            }
        }
        if ((*orb).steps % jump == 0) {
            writeInteData(*orb);
        }
    }
};

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
    crtbp() = default;
    ~crtbp() = default;
    std::ofstream summary;
    void inteSingleAdaptive(orbit3d &orbit, double endtime, size_t jump);
    void inteSingle(orbit3d &orbit, double endtime, size_t jump);
    void inteNbody(orbit3d orbits[], size_t n, double endtime, size_t jump);
    void writeSummary(orbit3d &orbit);
    static double jacobiConstant(const vec6 &x);
    static void eqOfMotion(const vec6 &x, vec6 &dxdt);
    static void eqOfVariation(const vec6 &x, vec6 &dxdt, const vec6 &p);
    static vec6 elementsToState(const vec6 &in);
    static vec6 inertialToRot(const vec6 &x, const double t);
    static vec6 rotToInertial(const vec6 &x, const double t);
    static vec6 elementsToRot(const vec6 &x, const double t);

  private:
    class crtbp_ode {
      public:
        void operator()(const vec6 &x, vec6 &dxdt, double t);
    };

    class crtbp_ode_variation {
      public:
        void operator()(const vec12 &x, vec12 &dxdt, double t);
    };

    static double true2mean(double theta, double e);
    static vec6 uxxMatrix(const vec3 &x);
};

static std::vector<std::pair<vec6, std::string>>
readInputFromTxt(const std::string &inputstring) {
    std::ifstream inputfile;
    inputfile.open(inputstring);
    std::vector<std::pair<vec6, std::string>> inputmatrix;
    std::string name;
    vec6 inputarray;
    if (inputfile.is_open()) {
        while (!inputfile.eof()) {
            for (size_t i = 0; i < 6; i++) {
                inputfile >> inputarray[i];
            }
            inputfile >> name;

            inputmatrix.push_back(std::make_pair(inputarray, name));
        }
    }
    return inputmatrix;
}

static vec3 readInfoFromTxt(const std::string &infostring) {
    std::ifstream infofile;
    infofile.open(infostring);
    vec3 infoarray;
    if (infofile.is_open()) {
        while (!infofile.eof()) {
            for (size_t i = 0; i < infoarray.size(); i++) {
                infofile >> infoarray[i];
            }
        }
    }
    return infoarray;
}

#endif