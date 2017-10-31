//
//  crtbp.cpp
//
//
//  Created by Protoss Probe on 2017/06/07.
//  Copyright © 2016-2017 probe. All rights reserved.
//

#ifndef _CRTBP_HPP_
#define _CRTBP_HPP_

#include "crtbp.hpp"
#include "utility.hpp"
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
const double pi2 = 2 * M_PI;
const double pi180 = pi / 180;
const double mu = 0.001;
const size_t MAX_NUMBER = 1000;
const double year = 365.25636042;
const std::string GLOBAL_OUTPUT_LOCATION = "assets/_output/";
const std::string GLOBAL_LOCATION = "assets/";
typedef boost::array<double, 2> vec2;
typedef boost::array<double, 3> vec3;
typedef boost::array<double, 4> vec4;
typedef boost::array<double, 5> vec5;
typedef boost::array<double, 6> vec6;
typedef boost::array<double, 7> vec7;
typedef boost::array<double, 12> vec12;

class orbit3d {
  private:
    double div = 1e-8;
    double megno_temp = 0;
    double megno_fast = 0;
    double megno_sum = 0;

  public:
    std::string name = "default";
    double jacobi0 = 0;
    double jacobi = 0;
    double time = 0;
    double dt = 0.001;
    double megno = 0;
    double megno_max = 0;
    size_t steps = 0;
    std::ofstream outputfile;
    vec12 vec = {{0, 0, 0, 0, 0, 0, div, div, div, div, div, div}};
    vec6 ele = {{1, 0, 0, 0, 0, 0}};
    vec6 vec_inertial = {{0, 0, 0, 0, 0, 0}};
    double getTime() { return time; }
    double getTimeYear() { return time / pi2; }
    vec6 getState() {
        return {{vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]}};
    }

    vec6 getDelta() {
        return {{vec[6], vec[7], vec[8], vec[9], vec[10], vec[11]}};
    }
    vec3 getPos() { return {{vec[0], vec[1], vec[2]}}; }
    void calInerState();
    vec3 getInerPos() {
        return {{vec_inertial[0], vec_inertial[1], vec_inertial[2]}};
    }
    vec3 getInerVel() {
        return {{vec_inertial[3], vec_inertial[4], vec_inertial[5]}};
    }
    vec6 getInerState() { return vec_inertial; }

    vec3 getVel() { return {{vec[3], vec[4], vec[5]}}; }
    std::string getName() { return name; };
    vec6 getElement() { return ele; };

    void setState(const vec6 &state);
    void setElement(vec6 element) { ele = element; };
    void setName(std::string newname) { name = newname; };
    void setDt(double t) { dt = t; };
    double errorJacobi();
    void setOutputFile() {
        std::string location = GLOBAL_OUTPUT_LOCATION + "Ast_" + getName();
        outputfile.open(location + ".txt");
    };
    void closeOutputFile() { outputfile.close(); };
    double getJacobi();
    vec2 deltaNormCal();

    double getLCN();
    void updateMEGNO(const double ticktime);
    double getMEGNO();
    double getMEGNOMax();
};

static void printInteData(orbit3d &orb) {
    std::cout << "State:  ";
    double j;
    BOOST_FOREACH (j, orb.getPos())
        std::cout << std::setw(12) << j << " ";
    std::cout << "  |   "
              << "Time: " << std::setw(6) << orb.getTime()
              << "   |   Jacobi: " << std::setw(8) << orb.getJacobi()
              << "   |   LCN: " << std::setw(8) << orb.getLCN()
              << "   |   MEGNO: " << std::setw(8) << orb.getMEGNO()
              << std::endl;
}

static void writeInteData(orbit3d &orb) {
    orb.calInerState();
    orb.outputfile << std::setw(8) << orb.getTimeYear() << '\t';
    double j;
    BOOST_FOREACH (j, orb.getInerState())
        orb.outputfile << std::setw(10) << j << '\t';
    orb.outputfile << std::setw(10) << orb.errorJacobi() << '\t'
                   << std::setw(10) << orb.getMEGNO() << std::endl;
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
            (*orb).getJacobi();
            auto megno = (*orb).getMEGNOMax();
            if (megno >= 8) {
                std::cout << megno << std::endl;
                throw 8.0;
                return;
            }
        }
        if ((*orb).steps % jump == 0) {
            writeInteData(*orb);
        }
    }
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