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

class CRTBP {
  public:
    CRTBP() = default;
    ~CRTBP() = default;
    CRTBP(double mu) : mu(mu) {}
    CRTBP(double mu, char option) : mu(mu), option(option) {}
    const double mu = 0.001;
    const char option = 'p';
    vec4 vector = {{0, 0, 0, 0}};
    vec4 elements = {{0, 0, 0, 0}};
    ofstream outputfile1, outputfile2, outputfile3, outputfile4, outputfile5;

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

    double jacobi_inte_ele(const vec4 &elements) {
        double a = elements[0];
        double e = elements[1];
        double result = 0.;
        if (option == 'p') {
            result = 1 / a + 2 * sqrt(a * (1 - e * e));
        } else if (option == 'r') {
            result = 1 / a - 2 * sqrt(a * (1 - e * e));
        }
        return result;
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
            stepper.do_step(planar_crtbp_ode(mu), x, 0.0, dt);
            current_t += dt;
            if (rot2r(x_temp) < 0 and rot2r(x) > 0) {
                write_keypoint(x, current_t);
            }
            // if (int(current_t / dt) % 50 == 0) {
            //     write_keypoint(x, current_t);
            // }
            write_detail(x, current_t);
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
        double phi = (elements[2] - elements[3] + t);
        double psi = atan2(x[1], x[0]);
        if (phi > pi or phi < -pi) {
            int n = floor(phi / (2 * pi));
            phi -= n * 2 * pi;
        }
        outputfile4 << elements[0] << '\t' << elements[1] << '\t' << phi << '\t'
                    << psi << endl;
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
        double E2 = jacobi_inte_ele(elements);
        outputfile1 << E << '\t' << endl;
        outputfile2 << endl;
        outputfile3 << E2 << '\t' << endl;

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

    double rot2r(vec4 x) { return x[0] * x[2] + x[1] * x[3]; }

    vec4 vector2elements(const vec4 &x) {
        // input: {{x,y,vx,vy}}
        // output: {{a,e,M,omega}}
        double r = sqrt(x[0] * x[0] + x[1] * x[1]);
        double h = x[0] * x[3] - x[2] * x[1];
        double angle = atan2(x[1], x[0]);
        double vr = x[2] * cos(angle) + x[3] * sin(angle);
        double vv = x[3] * cos(angle) - x[2] * sin(angle);
        double theta = atan2(vr, vv - (1 - mu) / h);
        double omega = angle - theta;
        // cout << angle << '\t' << theta << '\t' << omega << '\t' << endl;
        // cout << vr << '\t' << vv << endl;
        double e =
            sqrt(pow(vv * h / (1 - mu) - 1, 2) + pow(vr * h / (1 - mu), 2));
        double a = h * h / (1 - mu) / (1 - e * e);

        // double vr = (x[0] * x[2] + x[1] * x[3]) / r;
        // double v2 = x[2] * x[2] + x[3] * x[3];
        // double vr2 = vr * vr;

        // double co1 = v2 - (1 - mu) / r;
        // double co2 = r * vr;
        // double ex = 1 / (1 - mu) * (co1 * x[0] - co2 * x[2]);
        // double ey = 1 / (1 - mu) * (co1 * x[1] - co2 * x[3]);
        // double omega = atan2(ey, ex);
        // double e = sqrt(ex * ex + ey * ey);
        // double a = h * h / (1 - mu) / (1 - e * e);
        // double theta;
        // double temp = (h * h / ((1 - mu) * r) - 1) / e;
        // if (temp >= 1.) {
        //     theta = 0.0;
        // } else {
        //     theta = acos(temp);
        // }
        // if (vr < 0) {
        //     theta = 2 * pi - theta;
        // }
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
        vx = vr * cos(angle) - vv * sin(angle);
        vy = vr * sin(angle) + vv * cos(angle);
        if (option == 'p') {
            return {{x, y, vx, vy}};
        } else if (option == 'r') {
            return {{x, y, -vx, -vy}};
        } else {
            return {{0, 0, 0, 0}};
        }
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

    double F1_func(const vec4 &x, const double t) {
        double dist = sqrt(pow(x[0] - 1 + mu, 2) + pow(x[1], 2));
        auto y1 = rot2init(x, t);
        auto y2 = rot2init({{1 - mu, 0, 0, 0}}, t);
        double T1 = (y1[2] * y2[2] + y1[3] * y2[3]) / (1 - mu);
        return -1 / dist + T1;
    }

    double average(double a, double e, int limit) {
        // outputfile1.open("../assets/rot.txt");
        // outputfile2.open("../assets/init.txt");
        // outputfile3.open("../assets/elements.txt");
        elements = {{a, e, 0, 0}};
        cal_init_vec();
        double endtime = 12;
        double current_t = 0.0;
        double dt = 0.001;
        auto x = vector;
        vec4 x_temp;
        int cross = 0;
        double F1, result;

        runge_kutta_dopri5<vec4> stepper;
        while (current_t < endtime) {
            x_temp = x;
            F1 = F1_func(x, current_t);
            stepper.do_step(planar_crtbp_ode(mu), x, 0.0, dt);
            current_t += dt;
            result += F1 * dt;
            // write_detail(x, current_t);
            if (x_temp[1] * x[1] < 0) {
                cross++;
            }
            if (cross == limit)
                break;
        }

        // outputfile1.close();
        // outputfile2.close();
        // outputfile3.close();
        return abs(result) / (2 * pi);
    }

    void test() {
        double t = 0.0;
        print_vec(elements);
        cal_init_vec();
        auto x = vector;
        print_vec(x);
        auto x2 = rot2init(x, t);
        print_vec(x2);
        auto elements = vector2elements(x2);
        print_vec(elements);

        // cout << "Another Test !!!!" << endl;
        // vec4 test = {{1, 0, 0, 1}};
        // auto test2 = init2rot(test, 0);
        // test = rot2init(test2, 0);
        // print_vec(test);
        // print_vec(test2);
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

class DistFunc {
  public:
    DistFunc() = default;
    ~DistFunc() = default;
    DistFunc(int N) : N(N) {}
    const int N = 30;
    const int p = -1;
    const int q = 2;
    ofstream newcombfile;
    vector<double> newcomb_data;

    const vec31 ck = {
        {1.000000169188453,    -0.5000031057687094, 0.3749128690132108,
         -0.3119232338556135,  0.2808227070595242,  -0.2781248131852969,
         -0.01909432719811797, 0.6255846711707229,  4.391177347082226,
         -12.56056452779317,   -42.47138109688133,  115.3728242713690,
         277.3887782266406,    -722.4590131496257,  -1200.260700915786,
         3136.220460368545,    3516.623442481088,   -9652.445516335719,
         -6883.643191094829,   21198.48869866687,   8447.196762739993,
         -33009.87584926714,   -4999.308625893687,  35608.83781472251,
         -1613.928613838902,   -25316.80382412250,  5168.435382611849,
         10678.03836441569,    -3565.463660214670,  -2026.697767486133,
         894.4599340733893}};

    // for alpha = 1
    // const vec61 ak = {
    //     {1.,        -0.5,       0.375,     -0.3125,    0.273438,  -0.246094,
    //      0.225586,  -0.209473,  0.196381,  -0.185471,  0.176197,  -0.168188,
    //      0.16118,   -0.154981,  0.149446,  -0.144464,  0.13995,   -0.135834,
    //      0.132061,  -0.128585,  0.125371,  -0.122386,  0.119604,  -0.117004,
    //      0.114567,  -0.112275,  0.110116,  -0.108077,  0.106147,  -0.104317,
    //      0.102578,  -0.100924,  0.0993468, -0.0978415, 0.0964027, -0.0950255,
    //      0.0937057, -0.0924394, 0.0912231, -0.0900535, 0.0889279, -0.0878434,
    //      0.0867976, -0.0857884, 0.0848135, -0.0838711, 0.0829595, -0.0820769,
    //      0.081222,  -0.0803932, 0.0795892, -0.078809,  0.0780512, -0.0773148,
    //      0.076599,  -0.0759026, 0.0752249, -0.074565,  0.0739222, -0.0732958,
    //      0.072685}};

    // for alpha = 0.6299605249;
    const vec61 ak = {
        {5.60062058195,  -42.0386966482, 247.791328250,  -983.178765857,
         2578.54785797,  -4196.49169944, 3283.61033602,  942.542701585,
         -3347.44038646, -300.121174969, 3012.84399044,  1103.87714622,
         -2286.99900724, -2353.64980041, 476.876777623,  2586.48253552,
         1871.52323752,  -613.581961436, -2413.97244930, -2090.24479306,
         -128.226519777, 1841.96648376,  2435.43226398,  1382.48679593,
         -508.471706433, -2027.38964378, -2328.70975827, -1334.18433136,
         342.485470307,  1802.67377290,  2339.21224982,  1744.31554655,
         350.478893987,  -1180.91790935, -2170.23294409, -2212.49616640,
         -1316.53799019, 123.649043331,  1509.66253474,  2278.97176222,
         2130.88464696,  1122.85280657,  -345.942968151, -1693.32013823,
         -2375.98079317, -2095.63403629, -933.121215491, 665.181910900,
         2029.36315399,  2534.07305889,  1864.47558945,  226.123252164,
         -1652.54683551, -2750.16374043, -2244.07973656, -104.391912999,
         2464.16697942,  3191.69523532,  181.441911606,  -4519.08106038,
         1878.85824712}};

    double newcomb_operator_bk(int a, int b, int c, int d) {
        double result;
        if (c < 0 or d < 0) {
            result = 0;
        } else if (c == 0 and d == 0) {
            result = 1;
        } else if (c == 1 and d == 0) {
            result = b - a * 0.5;
        } else if (d == 0) {
            result = 2 * (2 * b - a) * newcomb_operator_bk(a, b + 1, c - 1, 0) +
                     (b - a) * newcomb_operator_bk(a, b + 2, c - 2, 0);
            result /= 4 * c;
        } else {
            result =
                -2 * (2 * b + a) * newcomb_operator_bk(a, b - 1, c, d - 1) -
                (b + a) * newcomb_operator_bk(a, b - 2, c, d - 2) -
                (c - 5 * d + 4 + 4 * b + a) *
                    newcomb_operator_bk(a, b, c - 1, d - 1);
            for (int j = 2; j <= min(c, d); j++) {
                result += 2 * (c - d + b) * pow(-1, j) * binomial(1.5, j) *
                          newcomb_operator_bk(a, b, c - j, d - j);
            }
            result /= 4 * d;
        }
        return result;
    }

    double newcomb_operator(int a, int b, int c, int d) {
        if (c < 0 or d < 0) {
            return 0;
        }
        int nn = 64;
        int nn1 = nn + 1;
        int at = a + nn;
        int bt = b + nn;
        return newcomb_data[d + c * nn1 + bt * nn1 * nn1 +
                            at * (2 * nn + 1) * nn1 * nn1];
    }

    void generate_newcomb() {
        int n = 10;
        // static double data[129][261][65][65];
        static double data[21][45][11][11];
        int at, bt;
        double result, temp1, temp2, temp3;
        temp1 = 0.0;
        temp2 = 0.0;
        temp3 = 0.0;
        for (int a = -n; a <= n; a++) {
            for (int b = -2 * n - 2; b <= 2 * n + 2; b++) {
                at = a + n;
                bt = b + 2 * n + 2;
                data[at][bt][0][0] = 1;
                data[at][bt][1][0] = b - a * 0.5;
            }
        }

        for (int c = 2; c <= n; c++) {
            for (int a = -n; a <= n; a++) {
                for (int b = -2 * n; b <= 2 * n; b++) {
                    at = a + n;
                    bt = b + 2 * n + 2;
                    data[at][bt][c][0] =
                        2 * (2 * b - a) * data[at][bt + 1][c - 1][0] +
                        (b - a) * data[at][bt + 2][c - 2][0];
                    data[at][bt][c][0] /= 4 * c;
                }
            }
        }
        for (int d = 1; d <= n; d++) {
            cout << "d: " << d << endl;
            for (int c = 0; c <= n; c++) {
                for (int a = -n; a <= n; a++) {
                    for (int b = -2 * n; b <= n; b++) {
                        at = a + n;
                        bt = b + 2 * n + 2;
                        temp1 = (d - 2 >= 0) ? data[at][bt - 2][c][d - 2] : 0;
                        temp3 = (c - 1 >= 0) ? data[at][bt][c - 1][d - 1] : 0;
                        data[at][bt][c][d] =
                            -2 * (2 * b + a) * data[at][bt - 1][c][d - 1] -
                            (b + a) * temp1 -
                            (c - 5 * d + 4 + 4 * b + a) * temp3;
                        for (int j = 2; j <= min(c, d); j++) {
                            temp2 = (c - j >= 0 and d - j >= 0)
                                        ? data[at][bt][c - j][d - j]
                                        : 0;
                            data[at][bt][c][d] += 2 * (c - d + b) * pow(-1, j) *
                                                  binomial(1.5, j) * temp2;
                        }
                        data[at][bt][c][d] /= 4 * d;
                    }
                }
            }
        }
        newcombfile.open("../assets/newcombfile_10.dat",
                         ios::out | ios::binary);
        for (int a = -n; a <= n; a++) {
            for (int b = -n; b <= n; b++) {
                for (int c = 0; c <= n; c++) {
                    for (int d = 0; d <= n; d++) {
                        at = a + n;
                        bt = b + 2 * n + 2;
                        // newcombfile << setprecision(14)
                        //             << (char *)&data[at][bt][c][d] << endl;
                        newcombfile.write((char *)&data[at][bt][c][d],
                                          sizeof(double));
                    }
                }
            }
        }
        newcombfile.close();
    }

    void read_newcomb() {
        string line;
        ifstream myfile;
        myfile.open("../assets/newcombfile_64.dat", ios::in);
        double data;
        while (!myfile.eof()) {
            myfile.read((char *)&data, sizeof(double));
            newcomb_data.push_back(data);
            if (myfile.fail()) {
                break;
            }
        }
        myfile.close();
    }

    void simple_test() {
        int n = 10;
        int nn = 64;
        int nn1 = nn + 1;
        int at, bt;
        for (int a = -n; a <= n; a++) {
            for (int b = -n; b <= n; b++) {
                for (int c = 0; c <= n; c++) {
                    for (int d = 0; d <= n; d++) {
                        at = a + nn;
                        bt = b + nn;
                        cout << setprecision(12)
                             << newcomb_operator(a, b, c, d) -
                                    newcomb_operator_bk(a, b, c, d)
                             << endl;
                    }
                }
            }
        }
    }

    double Y_coefficient(int n, int k, int s, int j) {
        int diff = j - k;
        int u1 = max(0, diff);
        int u2 = max(0, -diff);
        // cout << "Newcomb_Input: " << s + u1 << '\t' << s + u2
        // << endl;
        // cout << n << '\t' << k << '\t' << s + u1 << '\t' << s + u2 << '\t'
        //      << endl;
        return newcomb_operator(n, k, s + u1, s + u2);
    }

    double B_coefficient(int n, int k, int i, int m) {
        // int s_d1;
        // s_d1 = i - abs(k - m);
        // double s1;
        // if (abs(s_d1) % 2 == 1) {
        //     s1 = 0;
        // } else {
        //     s1 = Y_coefficient(n, k, s_d1 / 2, m);
        // }
        int s_d1, s_d2;
        s_d1 = i - abs(k - m);
        s_d2 = i - abs(k + m);
        double s1, s2;
        if (abs(s_d1) % 2 == 1) {
            s1 = 0;
        } else {
            s1 = Y_coefficient(n, k, s_d1 / 2, m);
        }
        if (abs(s_d2) % 2 == 1) {
            s2 = 0;
        } else {
            s2 = Y_coefficient(n, k, s_d2 / 2, -m);
        }
        if (m == 0) {
            return s1;
        }
        return s1 + s2;
    }

    double C_coefficient(int n, int k, int i, int m) {
        int s_d1, s_d2;
        s_d1 = i - abs(k - m);
        s_d2 = i - abs(k + m);
        double s1, s2;
        if (abs(s_d1) % 2 == 1) {
            s1 = 0;
        } else {
            s1 = Y_coefficient(n, k, s_d1 / 2, m);
        }
        if (abs(s_d2) % 2 == 1) {
            s2 = 0;
        } else {
            s2 = Y_coefficient(n, k, s_d2 / 2, -m);
        }
        if (m == 0) {
            return 0.0;
        }
        return s1 - s2;
    }

    double D_coefficient(int i, int j, int k, int m, int n, int l) {
        // i - alpha
        // j - e1
        // k - e2
        // l - delta_omega
        // m - M1
        // n - M2
        double co = 1 / (2. * gamma_m(m) * gamma_m(n));
        double result = co *
                        (B_coefficient(i, l, j, abs(m)) +
                         sign(m) * C_coefficient(i, l, j, abs(m))) *
                        (B_coefficient(-i - 1, l, k, abs(n)) +
                         sign(n) * C_coefficient(-i - 1, l, k, abs(n)));
        return result;
    }

    double A_coefficient(int k, int i) {
        double result = 0;
        double co = pow(-1, k) * gamma_m(k);
        for (int j = i; j <= min(2 * i, N - k); j++) {
            result += ck[j + k] * binomial(j + k, 2 * j - 2 * i + k) *
                      binomial(2 * j - 2 * i + k, j - i);
        }
        result *= co;
        return result;
    }

    double I_coefficient(int i, int j) {
        return 1 / (2. * gamma_m(j)) *
               (B_coefficient(1, 1, i, abs(j)) +
                sign(j) * C_coefficient(1, 1, i, abs(j)));
    }

    double R_coefficient(int i, int j, int k, int m, int n, int l) {
        double result = 0;
        if (abs(i - l) % 2 != 1 and (i - l) >= 0 and (i + l) <= 2 * N) {
            result = D_coefficient(i, j, k, m, n, l);
            if (abs(result) > 1e-12) {
                result *= A_coefficient(l, (i - l) / 2);
            }
        }
        if (l == 0) {
            result -= ak[i] * m * n * I_coefficient(j, m) * I_coefficient(k, n);
        }
        return result;
    }

    double R_bar_coefficient(int i, int j, int k, int u, int l) {
        double result;
        int m = u * p;
        int n = u * (p + q);
        return R_coefficient(i, j, k, m, n, l);
    }

    double binomial(double n, int m) {
        if (m == 0) {
            return 1;
        }
        double result = 1;
        for (int i = 0; i < m; i++) {
            result *= n - i;
        }
        result /= fact(m);
        return result;
    }

    double fact(int n) {
        double result;
        if (n == 1)
            return 1;
        result = fact(n - 1) * n;
        return result;
    }

    double gamma_m(int m) { return (m == 0) ? 0.5 : 1; }

    double true2mean(double theta, double e) {
        double E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(theta / 2));
        return E - e * sin(E);
    }

    void test() {
        // validate equation (23)
        double a = 1;
        double r, LHS, RHS, x;
        double e = 0.2;
        double M;
        for (double f = 0; f < 2 * pi; f += 0.1) {
            r = a * (1 - e * e) / (1 + e * cos(f));
            x = r * cos(f);
            LHS = x / a;
            RHS = 0;
            M = true2mean(f, e);
            for (int i = 0; i < 15; i++) {
                for (int j = -10; j <= 10; j++) {
                    RHS += pow(e, i) * I_coefficient(i, j) * cos(j * M);
                }
            }
            cout << setprecision(12) << "f: " << f << "\tLHS " << LHS
                 << "\terror: " << LHS - RHS << endl;
        }
    }

    void test2() {
        // validate equation (16) with some corrections
        int n = 5;
        int k = 2;
        double a = 1;
        double r, LHS, RHS;
        double e = 0.2;
        double M;
        for (double f = 0; f < 2 * pi; f += 0.1) {
            r = a * (1 - e * e) / (1 + e * cos(f));
            LHS = pow(r / a, n) * cos(k * f);
            RHS = 0;
            M = true2mean(f, e);
            for (int i = 0; i < 15; i++) {
                for (int m = 0; m < 10; m++) {
                    RHS += pow(e, i) * B_coefficient(n, k, i, m) * cos(m * M);
                }
            }
            cout << setprecision(12) << "f: " << f << "\tLHS " << LHS
                 << "\terror: " << LHS - RHS << endl;
        }
    }

    void test3() {
        // validate equation (34)
        int i = 2;
        int j = 2;
        int k = 2;
        int l = 0;
        double a_res = pow(0.5, 2. / 3.);
        double a = a_res;
        for (double e1 = 0.; e1 < 0.8; e1 += 0.01) {
            double result;
            double sum = 0;
            int m, n;
            double sigma1 = (90. / 180.) * pi;
            double sigma2 = 0.;
            for (int i = 0; i <= 5; i++) {
                for (int j = 0; j <= 4; j++) {
                    for (int k = 0; k <= 0; k++) {
                        for (int l = 0; l <= 4; l++) {
                            for (int u = -10; u <= 10; u++) {
                                result = R_bar_coefficient(i, j, k, u, l);
                                if (abs(result) > 1e-12) {
                                    m = u * 1;
                                    n = u * 2;
                                    sum += result * pow(a - a_res, i) *
                                           pow(e1, j) *
                                           cos((m - l) * sigma1 +
                                               (l - n) * sigma2);
                                    // cout << "i: " << i <<
                                    // '\t' << "j: "
                                    // << j
                                    //      << '\t' << "k: " <<
                                    //      k << '\t'
                                    //      << "l: " << l <<
                                    //      '\t' << "u: "
                                    //      << u
                                    //      << '\t' << "R: " <<
                                    //      result << endl;
                                }
                            }
                        }
                    }
                }
            }
            cout << "sum: " << sum << endl;
        }
    }

    void test4() {
        // validate equation (28)
        double mu = 0.001;
        CRTBP system(mu, 'p');
        double a = pow(0.5, 2. / 3.);
        // a = 1.0;
        double e1 = 0.4;
        double f = 0.0;
        double omega = 0;
        auto v1 = system.elements2vector({{1, 0, 0, 0}});
        auto v2 = system.elements2vector({{a, e1, f, omega}});

        double result = 0;
        double sum = 0;
        int num = 30;

        for (int j = 0; j <= 60; j++) {
            cout << j << endl;
            for (int k = 0; k <= 0; k++) {
                for (int m = -num; m <= num; m++) {
                    for (int n = -num; n <= num; n++) {
                        for (int i = 0; i <= 2 * N; i++) {
                            for (int l = 0; l <= N; l++) {
                                result = R_coefficient(i, j, k, m, n, l);
                                if (abs(result) > 1e-12) {
                                    sum += result * pow(a, i) * pow(e1, j);
                                }
                            }
                        }
                    }
                }
            }
        }
        double dist = sqrt(pow(v1[0] - v2[0], 2) + pow(v1[1] - v2[1], 2));
        result = 1 / dist - (v1[2] * v2[2] + v1[3] * v2[3]);
        cout << "error: " << sum - result << endl;
    }

    void test5() {
        // validate equation (17)
        double mu = 0.001;
        CRTBP system(mu, 'p');
        double a = pow(0.5, 2. / 3.);
        a = 0.6;
        double e1 = 0.2;
        double f = 0.3;
        double omega = 0;
        auto v1 = system.elements2vector({{1, 0, 0, 0}});
        auto v2 = system.elements2vector({{a, e1, f, omega}});

        double M = true2mean(f, e1);
        double result, sum, co1, co2;
        int q;
        sum = 0.;
        int num = 30;
        for (int l = 0; l <= N; l++) {
            for (int i = 0; i <= N - l; i++) {
                q = 2 * i + l;
                co1 = 2 * pow((1 - e1 * e1) / (1 + e1 * cos(f)), q) *
                      cos(l * (f + omega));
                co2 = 0;
                for (int j = 0; j <= 60; j++) {
                    for (int k = 0; k <= 0; k++) {
                        for (int m = -num; m <= num; m++) {
                            for (int n = -num; n <= num; n++) {
                                result = D_coefficient(q, j, k, m, n, l);
                                if (abs(result) > 1e-12) {
                                    co2 += result * pow(e1, j) *
                                           cos(m * M + l * omega);
                                }
                            }
                        }
                    }
                }
                cout << q << '\t' << l << '\t' << co1 - co2 << endl;
                // cout << "i: " << i << '\t' << "l: " << l <<
                // '\t' << "co: " << co
                //      << endl;
                sum += co2 * pow(a, q) * A_coefficient(l, i);
            }
        }
        double dist = sqrt(pow(v1[0] - v2[0], 2) + pow(v1[1] - v2[1], 2));
        cout << "error: " << sum - 1 / dist << endl;
    };

    void test6() {
        // validate equation (17)
        double a_res = pow(0.5, 2. / 3.);
        double a = a_res;
        double e1 = 0.4;
        double f = 0.;
        double M = true2mean(f, e1);
        double delta_omega = 0;
        double result, sum;

        int num = 30;
        int q = 4;
        int l = 1;
        for (int l = 0; l <= N; l++) {
            for (int i = 0; i <= N - l; i++) {
                q = 2 * i + l;
                sum = 0.;
                for (int j = 0; j <= 40; j++) {
                    for (int k = 0; k <= 0; k++) {
                        for (int m = -num; m <= num; m++) {
                            for (int n = -num; n <= num; n++) {
                                result = D_coefficient(q, j, k, m, n, l);
                                if (abs(result) > 1e-12) {
                                    sum += result * pow(e1, j) *
                                           cos(m * M + l * delta_omega);
                                }
                            }
                        }
                    }
                }
                cout << q << '\t' << l << '\t';
                cout << "error: "
                     << sum - 2 * pow((1 - e1 * e1) / (1 + e1 * cos(f)), q) *
                                  cos(l * (f + delta_omega))
                     << endl;
            }
        }
    };

    void test7() {
        // validate equation (13)
        double sum = 0.;
        double r1 = 0.5;
        int m;
        int N = 30;
        for (int k = 0; k <= N; k++) {
            for (int i = 0; i <= N - k; i++) {
                m = 2 * i + k;
                sum += pow(r1, m) * A_coefficient(k, i);
            }
        }
        cout << sum - 1 << endl;
    };

    void test8() {
        // validate equation (11)
        double sum = 0.;
        double co;
        double r1 = 0.5;
        double S = 0;
        double x;
        double rho = r1;
        x = rho * rho - 2 * rho * cos(S);
        for (int k = 0; k <= 30; k++) {
            double co = 0.;
            for (int j = 0; j <= k; j++) {
                co += pow(-2, j) * binomial(k, j) * pow(rho, 2 * k - j);
            }
            cout << "k: " << k << '\t' << "co: " << co << endl;
            cout << pow(x, k) << endl;
            // sum = ck[k] * k * pow(rho, 2 * k);
            // sum += ck[k] * pow(0.5, k);
            sum += ck[k] * pow(x, k);
        }
        cout << x << endl;
        cout << sum - 2 << endl;
    }

    void test10() {
        // validate equation (27)
        double a_res = pow(0.5, 2. / 3.);
        double a1 = a_res;
        double e1 = 0.5;
        double sum = 0;
        int N = 30;

        double co1, co2;
        for (int j = 0; j <= 15; j++) {
            for (int k = 0; k <= 0; k++) {
                for (int m = -20; m <= 20; m++) {
                    for (int n = -20; n <= 20; n++) {
                        co2 = 0;
                        for (int i = 0; i <= 2 * N; i++) {
                            co1 = 1 / sqrt(a1);
                            co2 += ak[i] * pow(a1, i);
                        }
                        sum += m * n * I_coefficient(j, m) *
                               I_coefficient(k, n) * pow(e1, j) * co2;
                    }
                }
            }
        }

        cout << "result: "
             << sum - 1 * (sqrt(1 / a1) * sqrt((1 + e1) / (1 - e1))) << endl;
    }
};

int main() {
    double mu = 0.001;
    CRTBP system(mu, 'p');
    double a_res = pow(0.5, 2. / 3.);
    // system.print_vec(system.elements2vector({{1, 0, 0, 0}}));
    // system.print_vec(system.elements2vector({{a_res, 0.2, 0,
    // 0}})); a_res = 1.0167; double e = 0.0; for (double e =
    // 0.0; e < 0.8; e += 0.01) {
    //     cout << system.average(a_res, e, 2) << ", ";
    // }

    // system.print_vec(system.vector);
    // // cout << "delta: " << system.deltax1x2(1.57, -0.45) <<
    // endl; system.integ_output(system.vector, 1000 * pi,
    // 0.005);
    // // system.test();
    // ofstream outputfile;
    // outputfile.open("../assets/B_C.txt");

    DistFunc func;
    // int k = 1;
    // int i = 1;
    // for (int k = 0; k <= 10; k++) {
    //     for (int i = -10; i <= 10; i++) {
    //         cout << "k: " << k << '\t' << "i: " << i << '\t';
    //         cout << "I: " << func.I_coefficient(k, i) <<
    //         endl;
    //     }
    // }
    // func.generate_newcomb();
    func.read_newcomb();
    // func.simple_test();
    func.test4();
    // double e = 0.1;
    // system.elements = {{a_res + 0.01, e, 0.0, 0}};
    // system.cal_init_vec();
    // cout << system.F1_func(system.vector, 0) << endl;
    // cout << "DONE!" << endl;
}

;
