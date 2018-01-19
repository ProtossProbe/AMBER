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

int main(int argc, char *argv[]) {
    string switchstring(argv[3]);
    ifstream switchfile;
    switchfile.open(switchstring);

    char swi;
    crtbp system;
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
        system.inteNbody(orbits, number, endt, jump);
        cout << "Time: " << omp_get_wtime() - start << endl;
        break;
    }

    case '2': {
        cout << "Execute Programm 2: Generate Single Averaged Hamiltonian"
             << endl;
        double N, Sz, S_max, S_min;
        size_t n;
        switchfile >> N;
        switchfile >> Sz;
        switchfile >> S_min;
        switchfile >> S_max;
        switchfile >> n;
        cout << "N constant is: " << N << endl << endl;
        cout << "Sz is: " << Sz << endl << endl;
        system.singleAverage(N, Sz, S_min, S_max, n);
        break;
    }
    case '3': {
        cout << "Execute Programm 3: Generate Double Averaged Hamiltonian"
             << endl;
        double H_val, a;
        switchfile >> H_val;
        switchfile >> a;
        cout << "H constant is: " << H_val << endl;
        cout << "a is: " << a << endl << endl;
        ;
        system.doubleAverage(H_val, a);
        break;
    }
    case '4': {
        cout << "Execute Programm 4: Find Planar Resonance Point" << endl;
        ofstream output;
        double Ns, S = 0.0058;
        output.open(GLOBAL_OUTPUT_LOCATION + "EqPoints.txt");
        for (Ns = -2.08; Ns <= -1.3; Ns += 0.001) {
            S = system.findEqPoint(Ns, 0, S);
            cout << setprecision(6) << "N constant is: " << Ns << endl << endl;
            cout << setprecision(9) << "S is: " << S << endl << endl;
            output << setprecision(12) << Ns << '\t' << S << endl;
        }
        output.close();
        break;
    }
    case '5': {
        cout << "Execute Programm 2: Generate Single Averaged omega Hamiltonian"
             << endl;
        double N, aa, Sz_max;
        size_t nn;
        switchfile >> N;
        switchfile >> aa;
        switchfile >> Sz_max;
        switchfile >> nn;
        cout << "N constant is: " << N << endl << endl;
        cout << "a is: " << aa << endl << endl;
        system.singleAverageOme(N, aa, Sz_max, nn);
        break;
    }
    }
    switchfile.close();
}