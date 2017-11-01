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
#include <exception>
#include <fstream>
#include <iostream>
#include <libiomp/omp.h>
#include <math.h>
#include <string>
#include <utility>
#include <vector>

// define some handy constants
const double pi = M_PI;
const double pi2 = 2 * M_PI;
const double pi180 = pi / 180;
const double pi_180 = 180 / pi;

/**
 * with the assumption of the total mass of two bodies is equal to 1,
 * we define:
 *      mu - mass of the planet (a.k.a the secondaty);
 *      mu_sun - mass of the center body (a.k.a the primary);
 *      year - days in one sidereal year (to convert unit between real model and
        dimensionless model);
 *      muu - 1 - mu;
 *      MAX_NUMBER - max number of particles can be handle in this program;
*/
const double mu = 0.001;
const double mu_sun = 1 - mu;
const double year = 365.25636042;
const double muu = mu_sun;
const size_t MAX_NUMBER = 1000;

// ouput file location
const std::string GLOBAL_OUTPUT_LOCATION = "assets/_output/";
const std::string GLOBAL_LOCATION = "assets/";

// define some handy vector types
typedef boost::array<double, 2> vec2;
typedef boost::array<double, 3> vec3;
typedef boost::array<double, 4> vec4;
typedef boost::array<double, 5> vec5;
typedef boost::array<double, 6> vec6;
typedef boost::array<double, 7> vec7;
typedef boost::array<double, 12> vec12;

// orbit3d class
class orbit3d {
  private:
    double div = 1e-8;
    double megno_temp = 0;
    double megno_fast = 0;
    double megno_sum = 0;
    double ticktime = 0.0;
    size_t new_steps = 0;

    vec2 deltaNormCal();
    void updateInerState();
    void updateJacobi();
    void updateMEGNO();
    void updateElements();

  public:
    // name of this orbit
    std::string name = "default";

    // initial jacobi constant
    double jacobi0 = 0;

    // current jacobi constant
    double jacobi = 0;

    // current relative jacobi error
    double jacobi_err = 0;

    // current time
    double time = 0;

    // current time interval
    double dt = 0.001;

    // current average megno (not yet calculated until time > ticktime)
    double megno = 0;

    // maximum megno value during the ingegration
    double megno_max = 0;

    // current number of steps;
    size_t steps = 0;

    // assigned file for output;
    std::ofstream outputfile;

    // a vector composed of a state vector (in the rotating frame) and a
    // variational vector
    vec12 vec = {{0, 0, 0, 0, 0, 0, div, div, div, div, div, div}};

    // an orbital elements vector {{a,e,I,g,n,f}} or {{a,e,I,g,n,m}}
    // a = semi-major axis (in AU)
    // e = eccentricity
    // I = inclination (degrees)
    // g = argument of pericentre (degrees)
    // n = longitude of the ascending node (degrees)
    // f = true anomaly (degrees)
    // m = mean anomaly (degrees)
    vec6 ele = {{1, 0, 0, 0, 0, 0}};

    // a state vector in the inertial frame
    vec6 vec_inertial = {{0, 0, 0, 0, 0, 0}};

    // get & set functions
    // declared and defined in this header file
    double getTime() { return time; }
    double getTimeYear() { return time / pi2; }
    vec6 getState() {
        return {{vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]}};
    }
    vec6 getDelta() {
        return {{vec[6], vec[7], vec[8], vec[9], vec[10], vec[11]}};
    }
    vec3 getPos() { return {{vec[0], vec[1], vec[2]}}; }
    vec3 getVel() { return {{vec[3], vec[4], vec[5]}}; }
    vec6 getInerState() { return vec_inertial; }
    vec3 getInerPos() {
        return {{vec_inertial[0], vec_inertial[1], vec_inertial[2]}};
    }
    vec3 getInerVel() {
        return {{vec_inertial[3], vec_inertial[4], vec_inertial[5]}};
    }
    std::string getName() { return name; };
    vec6 getElements() { return ele; };
    double getJacobi() { return jacobi; };
    double getJacobiError() { return jacobi_err; };
    double getMEGNO() { return megno; };
    double getMEGNOMax() { return megno_max; };
    void setState(const vec6 &state) {
        for (size_t i = 0; i < state.size(); i++) {
            vec[i] = state[i];
        }
    };
    void setElement(vec6 element) { ele = element; };
    void setName(std::string newname) { name = newname; };
    void setDt(double t) { dt = t; };
    void setOutputFile() {
        std::string location = GLOBAL_OUTPUT_LOCATION + "Ast_" + getName();
        outputfile.open(location + ".txt");
    };
    void setTickTime(double coef) { ticktime = coef * time; }
    void setInitial(vec6 elements, double dtt, std::string na);
    void closeOutputFile() { outputfile.close(); };

    // declared here, defined in crtbp.cpp
    double updatePerStep(const double t);
    void updatePerOutput();
    double getLCN();
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
    orb.outputfile << std::setw(8) << orb.getTimeYear() << '\t';
    double j;
    BOOST_FOREACH (j, orb.getElements())
        orb.outputfile << std::setw(10) << j << '\t';
    orb.outputfile << std::setw(10) << orb.getJacobiError() << '\t'
                   << std::setw(10) << orb.getMEGNO() << std::endl;
}

struct observer {
    orbit3d *orb;
    size_t jump;
    observer(orbit3d &orbit, size_t jump) : orb(&orbit), jump(jump) {}
    void operator()(const vec12 &x, double t) {
        if (t > 0) {
            auto megno = (*orb).updatePerStep(t);
            if (megno >= 8) {
                throw 8.0;
                return;
            }
        }
        if ((*orb).steps % jump == 0) {
            (*orb).updatePerOutput();
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
    static vec6 stateToElements(const vec6 &in, const char option = 'f');
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