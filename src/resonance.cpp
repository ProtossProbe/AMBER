//
//  crtbp.cpp
//
//
//  Created by Protoss Probe on 2017/06/07.
//  Copyright © 2016-2017 probe. All rights reserved.
//

#include "resonance.hpp"

using namespace std;
using namespace boost::math;
using namespace boost::numeric::odeint;
using namespace ProbeUtility;

void resonance::calAEI()
{
    double L;
    if (isPrograde)
    {
        L = k * (N - S - Sz) / ord;
        e = 1 - S / L;
        I = 1 - Sz / (L * e);
    }
    else
    {
        L = k * (S + Sz - N) / ord;
        e = 1 - S / L;
        I = Sz / (L * e) - 1;
    }
    e = sqrt(1 - e * e);
    I = acos(I) * pi_180;
    a = pow(L, 2);
}

void resonance::calAction()
{
    double L = sqrt(a);
    double temp_e = sqrt(1 - e * e);
    double G = L * temp_e;
    double H = G * cos(I * pi180);
    if (isPrograde)
    {
        N = L * (1 + ord / k) - H;
        S = L - G;
        Sz = G - H;
    }
    else
    {
        N = L * (1 - ord / k) + H;
        S = L - G;
        Sz = G + H;
    }
}

void resonance::calH0()
{
    double val;
    val = isPrograde ? N - S - Sz : S - N + Sz;
    H0 = -ord * ord / (2 * pow(k * val, 2)) - val * kj / ord;
}

void resonance::calPhi()
{
    double lambda;
    if (isPrograde)
    {
        lambda = M + ome + Ome;
        phi = kj * Mj - k * lambda - (kj - k) * (ome + Ome);
    }
    else
    {
        lambda = M + ome - Ome;
        phi = k * lambda - kj * Mj + ord * (Ome - ome);
    }
}

void resonance::calM()
{
    double lambda;
    if (isPrograde)
    {
        lambda = (kj * Mj - (kj - k) * (ome + Ome) - phi) / k;
        M = lambda - ome - Ome;
    }
    else
    {
        lambda = (phi + kj * Mj - ord * (Ome - ome)) / k;
        M = lambda - ome + Ome;
    }
}

void resonance::calTrue()
{
    f = crtbp::mean2true(fmod(M, 360) * pi180, e, 1e-10) * pi_180;
}

void resonance::displayInfo()
{

    cout << (isPrograde ? "Prograde!" : "Retrograde!") << endl;
    cout << "Ratio is: " << kj << "/" << k << endl;
    double a_c = pow(double(k) / kj, 2. / 3.);
    cout << "Accurate a is: " << a_c << endl;
    cout << "Accurate N is: "
         << (isPrograde ? (sqrt(a_c) * abs(k - kj) / k)
                        : (-sqrt(a_c) * double(k + kj) / k))
         << endl
         << endl;
}

double resonance::dotDisturb(const vec3 &v, const vec3 &r1)
{
    vec3 delta = {{v[0] - r1[0], v[1] - r1[1], v[2] - r1[2]}};
    double r1_norm = vecNorm(r1);
    double delta_norm = vecNorm(delta);
    if (delta_norm == 0 or r1_norm == 0)
    {
        return numeric_limits<double>::infinity();
    }
    else
    {
        return 1 / delta_norm - vecDot(v, r1) / pow(r1_norm, 3);
    }
}

double resonance::calDotDisturb()
{
    calTrue();
    vec6 vec = crtbp::elementsToRot({{a, e, I, ome, Ome, f}}, Mj * pi180);
    return dotDisturb({{vec[0], vec[1], vec[2]}}, {{1 - mu, 0, 0}});
}

double resonance::ringDisturb(const vec3 &v,
                              const double ap)
{
    double xr = sqrt(v[0] * v[0] + v[1] * v[1]), z2 = v[2] * v[2];
    double p2 = pow(xr + ap, 2) + z2, q2 = pow(xr - ap, 2) + z2;
    double p = sqrt(p2);
    double k = 1 - q2 / p2;
    if (k >= 1)
    {
        return numeric_limits<double>::infinity();
    }
    return 4 / p * comp_ellint_1(sqrt(k));
}

double resonance::calRingDisturb(string planets)
{
    calTrue();
    vec6 vec = crtbp::elementsToRot({{a, e, I, ome, Ome, f}}, 0);
    vec3 pos_p = {{vec[0], vec[1], vec[2]}};
    double result = 0;
    double mu_temp, d_temp;
    size_t ind;
    for (char &pl : planets)
    {
        ind = pl - '1';
        mu_temp = mu_group[ind];
        d_temp = d_group[ind];
        result += mu_temp * ringDisturb(pos_p, d_temp / PLANET_DISTANCE);
    }
    return result;
}

double resonance::inteOneCycle()
{
    double incr = 360.0 / num * k, disturb = 0;
    for (size_t i = 0; i < num; i++)
    {
        Mj = incr * i;
        calM();
        disturb += calDotDisturb();
    }
    // cout << phi << '\t' << a << '\t' << e << '\t' << Mj << '\t' << M << '\t'
    //      << setprecision(10) << disturb << endl;
    return mu * disturb / num;
}

double resonance::inteRing(string planets)
{
    double incr = 360.0 / num, disturb = 0;
    for (size_t i = 0; i < num; i++)
    {
        M = incr * i;
        disturb += calRingDisturb(planets);
    }
    // cout << a << '\t' << e << '\t' << Mj << '\t' << M << '\t'
    //      << setprecision(10) << disturb << endl;

    return mu * disturb / num;
}

double resonance::inteTwoCycle(double phiAmp)
{
    double disturb = 0;
    size_t numt = num / 4;
    for (size_t i = 0; i < numt; i++)
    {
        phi = sin((double)i / numt * pi2) * phiAmp + phi_c;
        disturb += inteOneCycle();
    }
    phi = phi_c;
    return disturb / num;
}

void resonance::averagePhi()
{
    ofstream output;
    output.open(GLOBAL_OUTPUT_LOCATION + "SingleAve_" + to_string(N) + "_" +
                to_string(Min) + "_" + to_string(Max) + ".txt");
    displayInfo();
    double dS = (Max - Min) / numY;
    double result;
    for (S = Min; S <= Max + 1e-10; S += dS)
    {
        cout << "S: " << S << endl;
        calAEI();
        calH0();
        cout << a << '\t' << e << '\t' << I << '\t' << endl;
        for (phi = -180; phi <= 180; phi += 1)
        {
            result = inteOneCycle() - H0;
            output << setprecision(12) << result << endl;
        }
    }
    output.close();
}

void resonance::averageOme(double phiAmp)
{
    double Sz_c;
    Sz_c = (isPrograde ? N - ord / k * sqrt(a) : N + ord / k * sqrt(a));
    Max = Sz_c;
    cout << "Critical Sz is: " << Sz_c << endl;

    ofstream output;
    output.open(GLOBAL_OUTPUT_LOCATION + "SingleAveOme_" + to_string(N) + "_" +
                to_string(Min) + "_" + to_string(Max) + ".txt");
    displayInfo();
    double result;
    double dSz = (Max - Min - 1e-8) / numY;
    for (Sz = Min; Sz <= Sz_c; Sz += dSz)
    {
        cout << "Sz: " << Sz << endl;
        S = (isPrograde ? N - ord / k * sqrt(a) - Sz
                        : N + ord / k * sqrt(a) - Sz);
        calAEI();
        calH0();
        cout << a << '\t' << e << '\t' << I << '\t' << phi << '\t' << endl;
        for (ome = -180; ome <= 180; ome += 2)
        {
            if (phiAmp == 0)
                result = inteOneCycle() - H0;
            else
            {
                result = inteTwoCycle(phiAmp) - H0;
            }
            output << setprecision(12) << result << endl;
        }
    }
    output.close();
}

void resonance::averageOmeRing(string planets)
{
    double cosi = cos(iMax * pi180);
    double eMax = sqrt(1 - cosi * cosi);
    double H = cosi;
    cout << "Maximum I is: " << iMax << endl;
    cout << "Maximum e is: " << eMax << endl;
    cout << "Planet You Choose: " << planets << endl;

    ofstream output;
    output.open(GLOBAL_OUTPUT_LOCATION + "RingAveOme_" + to_string(iMax) + "_" +
                to_string(a) + "_" + to_string(eMax) + ".txt");
    double result;
    double de = (eMax - Min - 1e-8) / numY;
    for (e = Min; e <= eMax; e += de)
    {
        I = acos(H / (sqrt(1 - e * e))) * pi_180;
        cout << "e: " << e << "\t";
        cout << "I: " << I << endl;
        for (ome = -180; ome <= 180; ome += 1)
        {
            result = inteRing(planets);
            output << setprecision(12) << result << endl;
        }
    }
    output.close();
}
