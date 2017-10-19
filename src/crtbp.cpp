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
#include <math.h>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace boost::math;
using namespace boost::numeric::odeint;

void crtbp::crtbp_ode::operator()(const vec6 &x, vec6 &dxdt, double t) {
    crtbp::eqOfMotion(x, dxdt);
}

void crtbp::crtbp_ode_variation::operator()(const vec12 &x, vec12 &dxdt,
                                            double t) {
    vec6 x1 = {{x[0], x[1], x[2], x[3], x[4], x[5]}};
    vec6 x2 = {{x[6], x[7], x[8], x[9], x[10], x[11]}};
    vec6 dxdt1, dxdt2;
    crtbp::eqOfMotion(x1, dxdt1);
    crtbp::eqOfVariation(x2, dxdt2, x1);
    dxdt = joinVector6(dxdt1, dxdt2);
};

void crtbp::eqOfMotion(const vec6 &x, vec6 &dxdt) {
    double x1 = x[0], x2 = x[1], x3 = x[2], x4 = x[3], x5 = x[4], x6 = x[5];
    double r1, r2;
    double xmu = x1 + mu, xmu1 = xmu - 1;
    double muu = 1 - mu;
    double x22 = x2 * x2, x33 = x3 * x3;
    r1 = sqrt(pow(xmu, 2) + x22 + x33);
    r2 = sqrt(pow(xmu1, 2) + x22 + x33);
    r1 = 1 / pow(r1, 3);
    r2 = 1 / pow(r2, 3);

    dxdt[0] = x4;
    dxdt[1] = x5;
    dxdt[2] = x6;
    dxdt[3] = -muu * xmu * r1 - mu * xmu1 * r2 + x1 + 2 * x5;
    dxdt[4] = -muu * x2 * r1 - mu * x2 * r2 + x2 - 2 * x4;
    dxdt[5] = -muu * x3 * r1 - mu * x3 * r2;
}

void crtbp::eqOfVariation(const vec6 &x, vec6 &dxdt, const vec6 &p) {
    double x1 = x[0], x2 = x[1], x3 = x[2], x4 = x[3], x5 = x[4], x6 = x[5];
    vec6 ux = crtbp::uxxMatrix({{p[0], p[1], p[2]}});
    double uxx = ux[0], uxy = ux[1], uxz = ux[2], uyy = ux[3], uyz = ux[4],
           uzz = ux[5];
    dxdt[0] = x4;
    dxdt[1] = x5;
    dxdt[2] = x6;
    dxdt[3] = x1 * uxx + x2 * uxy + x3 * uxz + 2 * x5;
    dxdt[4] = x1 * uxy + x2 * uyy + x3 * uyz - 2 * x4;
    dxdt[5] = x1 * uxz + x2 * uyz + x3 * uzz;
}

double crtbp::jacobiConstant(const vec6 &x) {
    double x1 = x[0], x2 = x[1], x3 = x[2], x4 = x[3], x5 = x[4], x6 = x[5];
    double xmu = x1 + mu, xmu1 = xmu - 1;
    double x22 = x2 * x2, x33 = x3 * x3;
    double r1, r2;
    r1 = sqrt(pow(xmu, 2) + x22 + x33);
    r2 = sqrt(pow(xmu1, 2) + x22 + x33);
    return x1 * x1 + x22 - x4 * x4 - x5 * x5 - x6 * x6 + 2 * (1 - mu) / r1 +
           2 * mu / r2;
}

void crtbp::inteSingle(orbit3d &orbit, double endtime, double dt) {
    runge_kutta_dopri5<vec12> stepper;
    orbit.dt = dt;
    while (orbit.time <= endtime) {
        orbit.updateMEGNO();
        writeInteData(orbit);
        stepper.do_step(crtbp_ode_variation(), orbit.vec, 0.0, orbit.dt);
        orbit.time += orbit.dt;
    }
}

vec6 crtbp::uxxMatrix(const vec3 &x) {
    double x1 = x[0], x2 = x[1], x3 = x[2];
    double uxx, uxy, uxz, uyy, uyz, uzz;

    double xmu = x1 + mu, xmu1 = xmu - 1, muu = 1 - mu;
    double xmuu = xmu * xmu, xmu11 = xmu1 * xmu1;
    double x22 = x2 * x2, x33 = x3 * x3;
    double r1, r13, r15, r2, r23, r25;
    r1 = sqrt(xmuu + x22 + x33);
    r2 = sqrt(xmu11 + x22 + x33);

    r13 = pow(r1, 3);
    r15 = r13 * pow(r1, 2);
    r23 = pow(r2, 3);
    r25 = r23 * pow(r2, 2);

    r13 = 1 / r13;
    r15 = 1 / r15;
    r23 = 1 / r23;
    r25 = 1 / r25;

    double mur25 = 3 * mu * r25, mur23 = mu * r23, muur15 = 3 * muu * r15,
           muur13 = muu * r13;

    uxx = 1 + mur25 * xmu11 - mur23 + muur15 * xmuu - muur13;
    uxy = mur25 * x2 * xmu1 + muur15 * x2 * xmu;
    uxz = mur25 * x3 * xmu1 + muur15 * x3 * xmu;
    uyy = 1 + mur25 * x22 - mur23 + muur15 * x22 - muur13;
    uyz = mur25 * x2 * x3 + muur15 * x2 * x3;
    uzz = mur25 * x33 - mur23 + muur15 * x33 - muur13;

    return {{uxx, uxy, uxz, uyy, uyz, uzz}};
}

const double orbit3d::getJacobi() {
    return crtbp::jacobiConstant(orbit3d::getState());
}
const double orbit3d::deltaNorm() { return vector6Norm(orbit3d::getDelta()); }
const double orbit3d::deltaDotNorm() {
    vec6 result;
    crtbp::eqOfVariation(orbit3d::getDelta(), result, orbit3d::getState());
    return vector6Norm(result);
}

double orbit3d::getLCN() {
    if (time > 0)
        return log(orbit3d::deltaNorm() / sqrt(6) / div) / time;
    else
        return 0.0;
}

void orbit3d::updateMEGNO() {
    double incr = orbit3d::deltaDotNorm() / orbit3d::deltaNorm() * dt * time;
    megno_temp += incr;
}

double orbit3d::getMEGNO() {
    if (time > 0) {
        megno = megno_temp / time * 2;
        return megno;
    }
    return 0.0;
}

void orbit3d::setState(const vec6 &state) {
    for (size_t i = 0; i < state.size(); i++) {
        vec[i] = state[i];
    }
}