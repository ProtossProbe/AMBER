#include <boost/array.hpp>
#include <boost/assign/std/vector.hpp>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char *argv[]) {
    double list[] = {1, 2, 3, 4};
    cout << sizeof(list) / sizeof(list[0]) << endl;
}