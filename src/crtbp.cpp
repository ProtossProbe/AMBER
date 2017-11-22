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

void crtbp::crtbp_two_body::operator()(const vec6 &x, vec6 &dxdt, double t) {
    crtbp::eqOfTwoBody(x, dxdt);
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

void crtbp::eqOfTwoBody(const vec6 &x, vec6 &dxdt) {
    double x1 = x[0], x2 = x[1], x3 = x[2], x4 = x[3], x5 = x[4], x6 = x[5];
    double r1;
    double x22 = x2 * x2, x33 = x3 * x3;
    r1 = sqrt(pow(x1, 2) + x22 + x33);
    r1 = 1 / pow(r1, 3);

    dxdt[0] = x4;
    dxdt[1] = x5;
    dxdt[2] = x6;
    dxdt[3] = -x1 * r1 + x1 + 2 * x5;
    dxdt[4] = -x2 * r1 + x2 - 2 * x4;
    dxdt[5] = -x3 * r1;
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
    // auto stepper = make_controlled(1e-12, 1e-12, orbit.dt,
    // stepper_const);
    orbit.setOutputFile();
    cout << "Start: " << orbit.getName() << endl;
    try {
        integrate_const(stepper, eq, orbit.vec, 0., endtime, orbit.dt,
                        observer(orbit, jump));
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
            << " " << orbit.getElements()[1] << " " << orbit.getElements()[2]
            << " " << orbit.getMEGNOMax() << endl;
}

void crtbp::inteSingle(orbit3d &orbit, double endtime, size_t jump) {
    crtbp_ode_variation eq;
    runge_kutta_dopri5<vec12> stepper;

    orbit.setOutputFile();
    cout << "Start: " << orbit.getName() << endl;
    size_t steps = 0;
    while (orbit.current_t <= endtime) {
        if (steps % jump == 0) {
            writeInteData(orbit);
        }
        stepper.do_step(eq, orbit.vec, 0.0, orbit.dt);
        orbit.current_t += orbit.dt;
        steps++;
    }
    orbit.closeOutputFile();
    cout << "End: " << orbit.getName() << endl;
    cout << endl;
}

void crtbp::inteNbody(orbit3d orbits[], size_t n, double endtime, size_t jump) {
    summary.open(GLOBAL_OUTPUT_LOCATION + "summary.out");
#pragma omp parallel for schedule(dynamic, 4) num_threads(6)
    for (size_t i = 0; i < n; i++) {
        crtbp::inteSingleAdaptive(orbits[i], endtime, jump);
    }

    summary.close();
}

vec6 crtbp::uxxMatrix(const vec3 &x) {
    double x1 = x[0], x2 = x[1], x3 = x[2];
    double uxx, uxy, uxz, uyy, uyz, uzz;

    double xmu = x1 + mu, xmu1 = xmu - 1;
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
    // output: {{x,y,z,u,v,w}} (inerital)
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
        out[3 + i] = sqrt(mu_sun / p) *
                     (-sin(f) * PVector[i] + (cos(f) + e) * QVector[i]);
    }
    return out;
}

vec6 crtbp::stateToElements(const vec6 &in, const char option) {
    // input: {{x,y,z,u,v,w}} (inerital)
    // output: {{a,e,I,g,n,f}} or {{a,e,I,g,n,m}}
    // a = semi-major axis (in AU)
    // e = eccentricity
    // I = inclination (degrees)
    // g = argument of pericentre (degrees)
    // n = longitude of the ascending node (degrees)
    // f = true anomaly (degrees)
    // m = mean anomaly (degrees)
    int i;
    double a, e, I, g, n, f;
    double temp1, temp2;
    vec3 R = {{in[0], in[1], in[2]}};
    vec3 V = {{in[3], in[4], in[5]}};
    double radius = vecNorm(R);
    double vel = vecNorm(V);
    double vr = vecDot(R, V) / radius;
    vec3 unitR, unitV, hvector, unith, evector, unite, nvector, unitN, temp;
    vecDevide(unitR, R, radius);
    vecDevide(unitV, V, vel);
    vec3Cross(hvector, R, V);
    double hnorm = vecNorm(hvector);
    vecDevide(unith, hvector, hnorm);
    temp1 = vel * vel - mu_sun / radius;
    temp2 = radius * vr;
    for (i = 0; i < 3; i++)
        evector[i] = (temp1 * R[i] - temp2 * V[i]) / mu_sun;
    e = vecNorm(evector);
    bool isCircle = (fabs(e) <= 1e-15);
    a = hnorm * hnorm / (mu_sun * (1 - e * e));
    vecDevide(unite, evector, e);
    I = acos(unith[2]);
    unitN = {{-unith[1], unith[0], 0}};
    if (vecNorm(unitN) == 0) {
        n = 0;
        if (isCircle) {
            g = 0;
            f = atan2(unitR[1] * unith[2], unitR[0]);
        } else {
            vec3Cross(temp, unite, unitR);
            g = atan2(unite[1] * unith[2], unite[0]);
            f = atan2(vecDot(unith, temp), vecDot(unite, unitR));
        }
    } else {
        vec3Cross(temp, unitN, unitR);
        n = atan2(unith[0], -unith[1]);
        f = atan2(vecDot(unith, temp), vecDot(unitN, unitR));
        if (isCircle) {
            g = 0;
        } else {
            vec3Cross(temp, unitN, unite);
            g = atan2(vecDot(unith, temp), vecDot(unite, unitN));
            f = f - g;
        }
    }
    if (g < 0) {
        g += 2 * pi;
    }
    if (n < 0) {
        n += 2 * pi;
    }
    I *= pi_180;
    g *= pi_180;
    n *= pi_180;

    if (option == 'f') {
        if (f < 0) {
            f += 2 * pi;
        }
        f *= pi_180;
        return {{a, e, I, g, n, f}};
    }
    if (option == 'm') {
        double m = crtbp::true2mean(f, e);
        if (m < 0) {
            m += 2 * pi;
        }
        m *= pi_180;
        return {{a, e, I, g, n, m}};
    } else {
        throw "option must be 'f' or 'm'!";
    }
}

vec6 crtbp::inertialToRot(const vec6 &x, const double t) {
    double c = cos(t), s = sin(t);
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
    double c = cos(t), s = sin(t);
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
    if (theta == 0.0) {
        return 0.0;
    }
    double E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(theta / 2));
    return E - e * sin(E);
}

double crtbp::mean2true(double M, double e) {
    if (M == 0.0) {
        return 0.0;
    }
    double result, E;
    if (M < pi and M > 0) {
        E = M + e / 2;
    } else {
        E = M - e / 2;
    };
    double incr = 1;
    while (fabs(incr) > 1e-14) {
        incr = keplerIteration(E, e, M);
        E += incr;
    }
    result = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
    return result;
}

double crtbp::keplerIteration(double E, double e, double M) {
    return -(E - e * sin(E) - M) / (1 - e * cos(E));
}

bool crtbp::isCross(const vec6 &vec_ref, const vec6 &vec,
                    const vec6 &vec_last) {
    double cross_last = vec_last[0] * vec_ref[1] - vec_last[1] * vec_ref[0];
    double cross_now = vec_ref[0] * vec[1] - vec_ref[1] * vec[0];
    return cross_last * cross_now > 0;
}
bool crtbp::isPeri(const double vr, const double vr_last) {
    return (vr_last < 0) and (vr >= 0);
}
double crtbp::radialVel(const vec6 &v) {
    vec3 pos = {{v[0], v[1], v[2]}};
    vec3 vel = {{v[3], v[4], v[5]}};
    return vecDot(pos, vel);
}

double crtbp::calSingleAve(const double N, const double S, const double sigma) {
    // 1/-1 resonance
    double B = S + N;
    double C = S - N;
    double a = pow(C / 2, 2);
    double e = sqrt(1 - pow(B / C, 2));

    // 2/1 resonance
    // double B = N - S;
    // double C = N - 2 * S;
    // double a = pow(B, 2);
    // double e = sqrt(1 - pow(C / B, 2));
    if (e < 0 or e >= 1) {
        return 0;
    }

    double omega = -2 * sigma;
    double T = pi2;
    size_t num = 720;
    double incr = pi2 / num, incrt = T / num, M = 0, f = 0, disturb = 0, t = 0;

    vec6 ele, vec;
    for (size_t i = 0; i < num; i++) {
        f = crtbp::mean2true(M, e) * pi_180;
        ele = {{a, e, 180, omega, 0, f}};
        vec = crtbp::elementsToRot(ele, t);
        disturb +=
            crtbp::disturbFunc({{vec[0], vec[1], vec[2]}}, {{1 - mu, 0, 0}});
        t += incrt;
        M += incr;
    }
    return mu * disturb / num;
}

void crtbp::singleAverage(const double N) {
    ofstream output;
    output.open(GLOBAL_OUTPUT_LOCATION + "SingleAve_" + to_string(N) + ".txt");

    double S_min = 0, S_max = 0.04, dS = 0.0002;
    double H0, result;
    for (double S = S_min; S <= S_max; S += dS) {
        cout << "S: " << S << endl;
        H0 = crtbp::genH0(N, S, 0);
        for (double sig = -180; sig <= 180; sig += 1) {
            result = crtbp::calSingleAve(N, S, sig);
            output << setprecision(10) << result << '\t' << H0 << endl;
        }
    }

    output.close();
}

double crtbp::calDoubleAve(const double H, const double a, const double e,
                           const double g) {
    double I = acos(H / sqrt(a * (1 - e * e))) * pi_180;
    double T = pi2;
    size_t num = 720;
    double incr = pi2 / num, incrt = T / num, M = 0, f = 0, disturb = 0, t = 0;
    vec6 ele, vec;
    for (size_t i = 0; i < num; i++) {
        f = crtbp::mean2true(M, e) * pi_180;
        ele = {{a, e, I, g, 0, f}};
        vec = crtbp::elementsToRot(ele, t);
        disturb += crtbp::ringDisturb({{vec[0], vec[1], vec[2]}}, 1);
        t += incrt;
        M += incr;
    }
    return mu * disturb / num / pi2;
}

void crtbp::doubleAverage(const double H, const double a) {
    ofstream output;
    output.open(GLOBAL_OUTPUT_LOCATION + "DoubleAve_" + to_string(H) + ".txt");
    int e_m = sqrt(1 - H * H / a) * 100;
    double e_min = 0, e_max = e_m / 100.0, de = 0.001;
    cout << "emax: " << e_max << endl;
    double result;
    for (double e = e_min; e <= e_max; e += de) {
        cout << "e: " << e << endl;
        for (double g = 0; g <= 360; g += 1) {
            result = calDoubleAve(H, a, e, g);
            output << setprecision(10) << result << endl;
        }
    }
    output.close();
}

double crtbp::ringDisturb(const vec3 &v, const double ap) {
    double xr = sqrt(v[0] * v[0] + v[1] * v[1]), z2 = v[2]*v[2];
    double p2 = pow(xr + ap, 2) + z2, q2 = pow(xr - ap, 2) + z2;
    double p = sqrt(p2);
    double k = 1 - q2 / p2;
    if (k >= 1) {
        return numeric_limits<double>::infinity();
    }
    return 4 / p * comp_ellint_1(k);
}

double crtbp::disturbFunc(const vec3 &v, const vec3 &r1) {
    vec3 delta = {{v[0] - r1[0], v[1] - r1[1], v[2] - r1[2]}};
    double r1_norm = vecNorm(r1);
    double delta_norm = vecNorm(delta);
    if (delta_norm == 0 or r1_norm == 0) {
        return numeric_limits<double>::infinity();
    } else {
        return 1 / delta_norm - vecDot(v, r1) / pow(r1_norm, 3);
    }
}

double crtbp::genH0(const double N, const double S, const double Sz) {
    // 1 / -1 resonance
    double val = S - N + Sz;
    return -2 / pow(val, 2) - val / 2;

    // double val = N - S - Sz;
    // return -1 / (2 * pow(val, 2)) - 2 * val;
}

void orbit3d::updateInerState() {
    vec_inertial =
        crtbp::rotToInertial(orbit3d::getState(), orbit3d::getTime());
}

void orbit3d::updateJacobi() {
    jacobi = crtbp::jacobiConstant(orbit3d::getState());
    double error = jacobi - jacobi0;
    if (error == 0) {
        jacobi_err = jacobi0;
    } else {
        jacobi_err = error / jacobi0;
    }
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
    if (current_t > 0)
        return log(vecNorm(orbit3d::getDelta()) / sqrt(6) / div) / current_t;
    else
        return 0.0;
}

void orbit3d::updateMEGNO() {
    vec2 delta = orbit3d::deltaNormCal();
    if (current_t > 0) {
        double incr = delta[1] / delta[0] * dt * current_t;
        megno_temp += incr;
        megno_fast = megno_temp / current_t * 2;
        if (current_t > ticktime) {
            new_steps++;
            megno_sum += megno_fast;
            megno = megno_sum / new_steps;
            if (megno > megno_max)
                megno_max = megno;
        }
    }
}

void orbit3d::updateElements() {
    ele = crtbp::stateToElements(vec_inertial, 'm');
}

// values need to be updated per step:
// steps, time, megno, megno_max
double orbit3d::updatePerStep(const double t) {
    steps++;
    current_t = t;
    orbit3d::updateMEGNO();
    return orbit3d::getMEGNOMax();
}

// values need to be updated per output:
// jacobi, jacobi_err, vec_inertial, ele
void orbit3d::updatePerOutput() {
    orbit3d::updateJacobi();
    orbit3d::updateInerState();
    orbit3d::updateElements();
}

// values need to be set in the first time:
// ele, dt, name, vec, jacobi0
void orbit3d::setInitial(vec6 elements, double dtt, double endt, string na) {
    ele = elements;
    dt = dtt;
    name = na;
    auto ele_true = orbit3d::getElements();
    ele_true[5] = crtbp::mean2true(ele_true[5] * pi180, ele_true[1]) * pi_180;
    orbit3d::setState(crtbp::elementsToRot(ele_true, 0));
    orbit3d::setTickTime(0.15 * endt);
    jacobi0 = crtbp::jacobiConstant(orbit3d::getState());
    ;
}
