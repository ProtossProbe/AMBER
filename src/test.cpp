#include "crtbp.hpp"

using namespace std;
using namespace ProbeUtility;

int main(int argc, char *argv[]) {
    vec6 test = {{-270, -180, -90, 45, 330, 450}};
    vec6 temp;
    vecMultiply(temp, test, pi180);
}