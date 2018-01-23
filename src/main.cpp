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

int main(int argc, char *argv[]) {
    string switchstring(argv[3]);
    ifstream switchfile;
    switchfile.open(switchstring);

    size_t k, kj;
    double N, Sz, Max, Min, a;
    int phi, phiAmp;
    size_t numY, num;
    char swi;
    crtbp inte_system;
    switchfile >> swi;
    switch (swi) {
    case '1': {
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
        for (size_t i = 0; i < number; i++) {
            orbits[i].setInitial(input[i], dt, endt, to_string(i + 1));
        }
        double start = omp_get_wtime();
        inte_system.inteNbody(orbits, number, endt, jump);
        cout << "Time: " << omp_get_wtime() - start << endl;
        break;
    }

    case '2': {
        cout << "Execute Programm 2: Generate Single Averaged Hamiltonian"
             << endl
             << endl;
        switchfile >> kj;
        switchfile >> k;
        switchfile >> N;
        switchfile >> Sz;
        switchfile >> Min;
        switchfile >> Max;
        switchfile >> numY;
        switchfile >> num;
        cout << "N constant is: " << N << endl;
        cout << "Sz is: " << Sz << endl << endl;
        resonance aver_system(k, kj, N, Sz, Min, Max, numY, num);
        aver_system.averagePhi();
        break;
    }
    case '3': {
        cout << "Execute Programm 3: Generate Single Averaged omega Hamiltonian"
             << endl;
        switchfile >> kj;
        switchfile >> k;
        switchfile >> N;
        switchfile >> a;
        switchfile >> phi;
        switchfile >> Max;
        switchfile >> numY;
        switchfile >> num;
        cout << "N constant is: " << N << endl;
        cout << "a is: " << a << endl << endl;
        resonance aver_system(k, kj, N, a, phi, Max, numY, num);
        aver_system.averageOme(0);
        break;
    }
    case '4': {
        cout << "Execute Programm 4: Generate Double Averaged omega Hamiltonian"
             << endl;
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
        cout << "a is: " << a << endl << endl;
        resonance aver_system(k, kj, N, a, phi, Max, numY, num);
        aver_system.averageOme(phiAmp);
        break;
    }
    case '5': {
        cout << "Execute Programm 5: Generate Double Averaged Hamiltonian"
             << endl;
        double H_val, a;
        switchfile >> H_val;
        switchfile >> a;
        cout << "H constant is: " << H_val << endl;
        cout << "a is: " << a << endl << endl;
        ;
        inte_system.doubleAverage(H_val, a);
        break;
    }
        // case '4': {
        //     cout << "Execute Programm 4: Find Planar Resonance Point" <<
        //     endl; ofstream output; double Ns, S = 0.0058;
        //     output.open(GLOBAL_OUTPUT_LOCATION + "EqPoints.txt");
        //     for (Ns = -2.08; Ns <= -1.3; Ns += 0.001) {
        //         S = inte_system.findEqPoint(Ns, 0, S);
        //         cout << setprecision(6) << "N constant is: " << Ns << endl <<
        //         endl; cout << setprecision(9) << "S is: " << S << endl <<
        //         endl; output << setprecision(12) << Ns << '\t' << S << endl;
        //     }
        //     output.close();
        //     break;
        // }
    }
    switchfile.close();
}