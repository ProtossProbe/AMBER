
#include <cmath>
#include <iostream>
int main() {
    double hpi = std::acos(-1) / 2;
    std::cout << "K(0) = " << std::comp_ellint_1(0) << '\n'
              << "π/2 = " << hpi << '\n'
              << "K(0.5) = " << std::comp_ellint_1(sqrt(2) / 2) << '\n'
              << "F(0.5, π/2) = " << std::ellint_1(0.5, hpi) << '\n';
}