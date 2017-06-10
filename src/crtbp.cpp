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

using namespace std;
using namespace boost::math;
using namespace boost::numeric::odeint;

class CRTBP {
  public:
    CRTBP() = default;
    ~CRTBP() = default;
    CRTBP(double mu) : mu(mu) {}
    const double mu = 0.001;
    vec4 vector = {{0, 0, 0, 0}};
    vec4 elements = {{0, 0, 0, 0}};
    ofstream outputfile1, outputfile2, outputfile3, outputfile4;

    class planar_crtbp_ode {
        double mu;

      public:
        planar_crtbp_ode(double mu) : mu(mu) {}
        void operator()(const vec4 &x, vec4 &dxdt, double t) {
            double r1, r2;
            double xmu = x[0] + mu;
            r1 = sqrt(pow(xmu, 2) + pow(x[1], 2));
            r2 = sqrt(pow(xmu - 1, 2) + pow(x[1], 2));
            r1 = pow(r1, 3);
            r2 = pow(r2, 3);

            dxdt[0] = x[2];
            dxdt[1] = x[3];
            dxdt[2] =
                -(1 - mu) * xmu / r1 - mu * (xmu - 1) / r2 + x[0] + 2 * x[3];
            dxdt[3] = -(1 - mu) * x[1] / r1 - mu * x[1] / r2 + x[1] - 2 * x[2];
        }
    };

    double jacobi_inte(const vec4 &x) {
        double xmu = x[0] + mu;
        double y2 = pow(x[1], 2);
        double r1, r2;
        r1 = sqrt(pow(xmu, 2) + y2);
        r2 = sqrt(pow(xmu - 1, 2) + y2);
        return x[0] * x[0] + y2 - x[2] * x[2] - x[3] * x[3] +
               2 * (1 - mu) / r1 + 2 * mu / r2;
    }

    void integ_output(vec4 x, const double endtime, const double dt = 0.01) {
        outputfile1.open("../assets/rot.txt");
        outputfile2.open("../assets/init.txt");
        outputfile3.open("../assets/elements.txt");
        outputfile4.open("../assets/key.txt");

        runge_kutta_dopri5<vec4> stepper;
        double current_t = 0.0;
        vec4 x_temp;
        write_keypoint(x, current_t);
        while (current_t < endtime) {
            x_temp = x;
            // write_detail(x, current_t);
            stepper.do_step(planar_crtbp_ode(mu), x, 0.0, dt);
            // if (rot2r(x_temp) < 0 and rot2r(x) > 0)
            current_t += dt;
            if ((int)(current_t / dt) % 250 == 0)
                write_keypoint(x, current_t);
        }

        outputfile1.close();
        outputfile2.close();
        outputfile3.close();
        outputfile4.close();
    }

    void write_keypoint(const vec4 &x, const double t) {
        outputfile4 << t << '\t';
        auto y = rot2init(x, t);
        auto elements = vector2elements(y);
        // double lambda = elements[2] + elements[3];
        double lambda_p = t;
        double phi = 2 * (lambda_p - elements[3]) - elements[2];
        if (phi > pi or phi < -pi) {
            int n = floor(phi / (2 * pi));
            phi -= n * 2 * pi;
        }
        outputfile4 << elements[0] << '\t' << elements[1] << '\t' << phi
                    << endl;
    }

    void write_detail(const vec4 &x, const double t) {
        outputfile1 << t << '\t';
        outputfile2 << t << '\t';
        outputfile3 << t << '\t';
        auto y = rot2init(x, t);
        auto elements = vector2elements(y);
        for (auto ele : x) {
            outputfile1 << ele << '\t';
        };
        for (auto ele : y) {
            outputfile2 << ele << '\t';
        };
        for (auto ele : elements) {
            outputfile3 << ele << '\t';
        };
        double E = jacobi_inte(x);
        outputfile1 << E << '\t' << endl;
        outputfile2 << endl;
        outputfile3 << endl;

        // cout << '\t' << t << '\t' << x[0] << '\t' << x[1] << '\t' << E <<
        // '\t'
        //      << endl;
    }

    vec4 rot2init(vec4 x, const double t) {
        double c = cos(t);
        double s = sin(t);
        x[0] += mu;

        vec4 result;
        result[0] = c * x[0] - s * x[1];
        result[1] = s * x[0] + c * x[1];
        result[2] = -s * (x[0] + x[3]) + c * (x[2] - x[1]);
        result[3] = s * (x[2] - x[1]) + c * (x[0] + x[3]);

        return result;
    }

    vec4 init2rot(vec4 x, const double t) {
        double c = cos(t);
        double s = sin(t);
        vec4 result;
        result[0] = c * x[0] + s * x[1];
        result[1] = -s * x[0] + c * x[1];
        result[2] = -s * (x[0] - x[3]) + c * (x[2] + x[1]);
        result[3] = -s * (x[2] + x[1]) - c * (x[0] - x[3]);

        result[0] -= mu;
        return result;
    }

    vec4 vector2elements(const vec4 &x) {
        // input: {{x,y,vx,vy}}
        // output: {{a,e,M,omega}}
        double r = sqrt(x[0] * x[0] + x[1] * x[1]);
        double vr = (x[0] * x[2] + x[1] * x[3]) / r;
        double v2 = x[2] * x[2] + x[3] * x[3];
        double vr2 = vr * vr;
        double h = x[0] * x[3] - x[2] * x[1];
        double co1 = v2 - (1 - mu) / r;
        double co2 = r * vr;
        double ex = 1 / (1 - mu) * (co1 * x[0] - co2 * x[2]);
        double ey = 1 / (1 - mu) * (co1 * x[1] - co2 * x[3]);
        double omega = atan2(ey, ex);
        double e = sqrt(ex * ex + ey * ey);
        double a = h * h / (1 - mu) / (1 - e * e);
        double theta;
        double temp = (h * h / ((1 - mu) * r) - 1) / e;
        if (temp >= 1.) {
            theta = 0.0;
        } else {
            theta = acos(temp);
        }
        if (vr < 0) {
            theta = 2 * pi - theta;
        }
        double M = true2mean(theta, e);
        return {{a, e, M, omega}};
    }

    double true2mean(double theta, double e) {
        double E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(theta / 2));
        return E - e * sin(E);
    }

    vec4 elements2vector(const vec4 &elements) {
        // input: {{a,e,theta,omega}}
        // output:: {{x,y,vx,vy}}
        double a, e, theta, omega, h, r, angle, x, y, vx, vy, vv, vr;
        a = elements[0];
        e = elements[1];
        theta = elements[2];
        omega = elements[3];
        h = sqrt(a * (1 - mu) * (1 - e * e));
        r = a * (1 - e * e) / (1 + e * cos(theta));
        angle = omega + theta;
        x = r * cos(angle);
        y = r * sin(angle);
        vv = (1 - mu) / h * (1 + e * cos(theta));
        vr = (1 - mu) / h * e * sin(theta);
        vx = vr * cos(angle) + vv * sin(angle);
        vy = -vr * sin(angle) + vv * cos(angle);

        return {{x, y, vx, vy}};
    }

    void print_vec(const vec4 &x) {
        for (auto ele : x) {
            cout << ele << '\t';
        }
        cout << endl;
    }

    void cal_init_vec() {
        auto x = elements2vector(elements);
        vector = init2rot(x, 0);
    }

    void test() {
        double t = 0.0;
        print_vec(elements);
        auto x = elements2vector(elements);
        print_vec(x);
        auto test_vec = init2rot(x, t);
        print_vec(test_vec);
        auto elements = vector2elements(rot2init(test_vec, t));
        print_vec(elements);
    }

    double delta21(double a) {
        double fs1 = 0.38762717;
        double fd = -1.19049274;
        double mum = mu / (1 - mu);
        double n = sqrt((1 - mu) / pow(a, 3));
        double omega_dot = 2 * mum * n * a * fs1;
        double e_beta = 1.5 * pow(n, 3) * mum * mum * pow(a * fd, 2);
        double alpha = 2 * (n - 2 + omega_dot);
        return alpha * pow(4 / e_beta, 1. / 3.);
    }

    double Hamiltion21(double delta, double Phi, double phi) {
        return delta * Phi + Phi * Phi - 2 * sqrt(2 * Phi) * cos(phi);
    }

    double deltax1x2(double x1, double x2) {
        return -(0.5 * (x2 * x2 + x1 * x1) + 4 / (x1 + x2));
    }

  private:
};

int main() {
    double mu = 0.001;
    CRTBP system(mu);
    double a_res = pow(0.5, 2. / 3.);
    double e = 0.01;

    system.elements = {{a_res, e, 0, 0}};
    system.cal_init_vec();
    // cout << "delta: " << system.delta21(a_res) << endl;
    // cout << "delta: " << system.deltax1x2(2.21, 2.35) << endl;
    system.integ_output(system.vector, 500 * pi, 0.001);
    cout << "DONE!" << endl;
}
