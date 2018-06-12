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

int main(int argc, char *argv[])
{
    string switchstring(argv[3]);
    ifstream switchfile;
    switchfile.open(switchstring);

    size_t k, kj;
    bool isPrograde;
    double N, Sz, Max, Min, a, iMax;
    int phi, phiAmp;
    size_t numY, num;
    char swi;
    string planets;
    crtbp inte_system;
    switchfile >> swi;
    switch (swi)
    {
    case '1':
    {
        cout << "Execute Programm 1: Numerical Integration" << endl;
        string inputstring(argv[1]);
        string infostring(argv[2]);
        auto input = readInputFromTxt(inputstring);
        auto info = readInfoFromTxt(infostring);
        size_t number = input.size();
        cout << "Total Number: " << number << endl;

        double endt = info[0] * pi2, dt = info[1] / year * pi2;
        int jump = info[2];

        orbit3d *orbits = new orbit3d[MAX_NUMBER];
        for (size_t i = 0; i < number; i++)
        {
            orbits[i].setInitial(input[i], dt, endt, to_string(i + 1));
        }
        double start = omp_get_wtime();
        inte_system.inteNbody(orbits, number, endt, jump);
        cout << "Time: " << omp_get_wtime() - start << endl;
        break;
    }

    case '2':
    {
        cout << "Execute Programm 2: Generate Single Averaged Hamiltonian"
             << endl
             << endl;
        switchfile >> isPrograde;
        switchfile >> kj;
        switchfile >> k;
        switchfile >> N;
        switchfile >> Sz;
        switchfile >> Min;
        switchfile >> Max;
        switchfile >> numY;
        switchfile >> num;
        cout << "N constant is: " << N << endl;
        cout << "Sz is: " << Sz << endl
             << endl;
        resonance aver_system(isPrograde, k, kj, N, Sz, Min, Max, numY, num);
        aver_system.averagePhi();
        break;
    }
    case '3':
    {
        cout << "Execute Programm 3: Generate omega Hamiltonian inside MMR"
             << endl;
        switchfile >> isPrograde;
        switchfile >> kj;
        switchfile >> k;
        switchfile >> N;
        switchfile >> a;
        switchfile >> phi;
        switchfile >> phiAmp;
        switchfile >> Max;
        switchfile >> numY;
        switchfile >> num;
        cout << "N constant is: " << N << endl;
        cout << "a is: " << a << endl
             << endl;
        resonance aver_system(isPrograde, k, kj, N, a, phi, Max, numY, num);
        aver_system.averageOme(phiAmp);
        break;
    }
    case '4':
    {
        cout << "Execute Programm 4: Generate omega Hamiltonian outside MMR"
             << endl;
        switchfile >> iMax;
        switchfile >> a;
        switchfile >> planets;
        switchfile >> numY;
        switchfile >> num;
        resonance aver_system(iMax, a, numY, num);
        aver_system.averageOmeRing(planets);
        break;
    }
        switchfile.close();
    }
}