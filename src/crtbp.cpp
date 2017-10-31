//
//  crtbp.cpp
//
//
//  Created by Protoss Probe on 2017/06/07.
//  Copyright © 2016-2017 probe. All rights reserved.
//

#include "crtbp.hpp"

using namespace std;
using namespace boost::math;
using namespace boost::numeric::odeint;
using namespace ProbeUtility;

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

    joinVec(dxdt, dxdt1, dxdt2);
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

void crtbp::inteSingleAdaptive(orbit3d &orbit, double endtime, size_t jump) {
    crtbp_ode_variation eq;
    bulirsch_stoer<vec12> stepper(1e-12, 1e-12, orbit.dt, orbit.dt);
    // runge_kutta_dopri5<vec12> stepper_const;
    // auto stepper = make_controlled(1e-12, 1e-12, orbit.dt, stepper_const);
    double tick = endtime * 0.15;
    orbit.setOutputFile();
    cout << "Start: " << orbit.getName() << endl;
    try {
        integrate_const(stepper, eq, orbit.vec, 0., endtime, orbit.dt,
                        observer(orbit, jump, tick));
    } catch (double megno) {
        cout << "Stop when MEGNO is larger than 8" << endl;
        orbit.megno_max = megno;
    }
    writeSummary(orbit);
    orbit.closeOutputFile();
    cout << "End: " << orbit.getName() << endl;
    cout << endl;
}

void crtbp::writeSummary(orbit3d &orbit) {
    summary << orbit.getName() << " " << setprecision(5) << orbit.getJacobi()
            << " " << orbit.getElement()[1] << " " << orbit.getElement()[2]
            << " " << orbit.getMEGNOMax() << endl;
}

void crtbp::inteSingle(orbit3d &orbit, double endtime, size_t jump) {
    crtbp_ode_variation eq;
    runge_kutta_dopri5<vec12> stepper;
    auto stepper_controlled = make_controlled(1e-12, 1e-12, orbit.dt, stepper);

    orbit.setOutputFile();
    cout << "Start: " << orbit.getName() << endl;
    size_t steps = 0;
    double tick = endtime * 0.15;
    while (orbit.time <= endtime) {
        orbit.updateMEGNO(tick);
        if (steps % jump == 0) {
            writeInteData(orbit);
        }
        stepper.do_step(eq, orbit.vec, 0.0, orbit.dt);
        orbit.time += orbit.dt;
        steps++;
    }
    orbit.closeOutputFile();
    cout << "End: " << orbit.getName() << endl;
    cout << endl;
}

void crtbp::inteNbody(orbit3d orbits[], size_t n, double endtime, size_t jump) {
    summary.open(GLOBAL_LOCATION + "summary.out");
#pragma omp parallel for schedule(dynamic, 4) num_threads(6)
    for (size_t i = 0; i < n; i++) {
        crtbp::inteSingleAdaptive(orbits[i], endtime, jump);
    }

    summary.close();
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

vec6 crtbp::elementsToState(const vec6 &in) {
    // input: {{a,e,I,g,n,f}}
    // output:: {{x,y,z,u,v,w}} (inerital)
    // a = semi-major axis (in AU)
    // e = eccentricity
    // I = inclination (degrees)
    // g = argument of pericentre (degrees)
    // n = longitude of the ascending node (degrees)
    // f = true anomaly (degrees)
    double a = in[0], e = in[1], I = in[2] * pi180, g = in[3] * pi180,
           n = in[4] * pi180, f = in[5] * pi180;
    vec6 out;
    double p = a * fabs(1.0 - e * e);

    double sini, cosi, sing, cosg, sinn, cosn;
    sini = sin(I);
    cosi = cos(I);
    sing = sin(n);
    cosg = cos(n);
    sinn = sin(g);
    cosn = cos(g);

    vec3 HVector = {{sini * sing, -sini * cosg, cosi}};

    vec3 PVector = {{cosg * cosn - sing * sinn * cosi,
                     sing * cosn + cosg * sinn * cosi, sinn * sini}};

    // QVector=[-cosg*sinn-sing*cosn*cosi;-sing*sinn+cosg*cosn*cosi;cosn*sini];
    vec3 QVector;
    vec3Cross(QVector, HVector, PVector);
    double r = 0.0;
    r = p / (1.0 + e * cos(f));
    for (int i = 0; i < 3; i++) {
        out[i] = r * (cos(f) * PVector[i] + sin(f) * QVector[i]);
        out[3 + i] =
            sqrt(1 / p) * (-sin(f) * PVector[i] + (cos(f) + e) * QVector[i]);
    }
    return out;
}

vec6 crtbp::inertialToRot(const vec6 &x, const double t) {
    double c = cos(t);
    double s = sin(t);
    vec6 result;
    result[0] = c * x[0] + s * x[1];
    result[1] = -s * x[0] + c * x[1];
    result[2] = x[2];
    result[3] = -s * (x[0] - x[4]) + c * (x[3] + x[1]);
    result[4] = -s * (x[3] + x[1]) - c * (x[0] - x[4]);
    result[5] = x[5];
    result[0] -= mu;
    return result;
}

vec6 crtbp::rotToInertial(const vec6 &x, const double t) {
    double c = cos(t);
    double s = sin(t);
    double x0 = x[0] + mu;
    vec6 result;
    result[0] = c * x0 - s * x[1];
    result[1] = s * x0 + c * x[1];
    result[2] = x[2];
    result[3] = -s * (x0 + x[4]) + c * (x[3] - x[1]);
    result[4] = s * (x[3] - x[1]) + c * (x0 + x[4]);
    result[5] = x[5];
    return result;
}

vec6 crtbp::elementsToRot(const vec6 &x, const double t) {
    return crtbp::inertialToRot(crtbp::elementsToState(x), t);
}

double crtbp::true2mean(double theta, double e) {
    double E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(theta / 2));
    return E - e * sin(E);
}

void orbit3d::calInerState(){
    vec_inertial = crtbp::rotToInertial(orbit3d::getState(),orbit3d::getTime());
}

double orbit3d::getJacobi() {
    jacobi = crtbp::jacobiConstant(orbit3d::getState());
    return jacobi;
}

double orbit3d::errorJacobi() {
    if (jacobi0 != 0) {
        return (jacobi - jacobi0) / jacobi0;
    }
    return 0.0;
}

vec2 orbit3d::deltaNormCal() {
    vec2 delta;
    vec6 delta_vec = orbit3d::getDelta();
    delta[0] = vecNorm(delta_vec);
    vec6 deltadot_vec;
    crtbp::eqOfVariation(delta_vec, deltadot_vec, orbit3d::getState());
    delta[1] = vecDot(deltadot_vec, delta_vec) / delta[0];
    return delta;
}

double orbit3d::getLCN() {
    if (time > 0)
        return log(vecNorm(orbit3d::getDelta()) / sqrt(6) / div) / time;
    else
        return 0.0;
}

void orbit3d::updateMEGNO(double ticktime) {
    vec2 delta = orbit3d::deltaNormCal();
    if (time > 0) {
        double incr = delta[1] / delta[0] * dt * time;
        megno_temp += incr;
        megno_fast = megno_temp / time * 2;
        if (time > ticktime) {
            size_t new_steps = steps - ticktime / dt;
            megno_sum += megno_fast;
            megno = megno_sum / (steps - new_steps);
            if (megno > megno_max)
                megno_max = megno;
        }
    }
}

double orbit3d::getMEGNO() { return megno; }

double orbit3d::getMEGNOMax() { return megno_max; }

void orbit3d::setState(const vec6 &state) {
    for (size_t i = 0; i < state.size(); i++) {
        vec[i] = state[i];
    }
}